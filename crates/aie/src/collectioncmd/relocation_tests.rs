// Included in the collection test module to reuse its scientific fixtures.
mod relocation {
    use super::*;

    #[test]
    fn legacy_copies_are_verified_by_content_on_preflight_and_actual_open() {
        use std::io::Seek;
        let rooted = scientific_archive_fixture();
        let legacy = test_path("relocation-legacy", "aie");
        let copy = test_path("relocation-legacy-copy", "aie");
        let mut reader = evidence_io::format::SectionReader::open(&rooted).unwrap();
        let mut file = std::fs::File::create(&legacy).unwrap();
        file.write_all(evidence_io::format::MAGIC).unwrap();
        file.write_all(&evidence_io::format::SEEKABLE_VERSION.to_le_bytes())
            .unwrap();
        let mut directory = Vec::new();
        for (name, _, raw_len, _) in reader.entries().to_vec() {
            let (compressed, _) = reader.read_compressed(&name).unwrap();
            let offset = file.stream_position().unwrap();
            file.write_all(&[name.len() as u8]).unwrap();
            file.write_all(name.as_bytes()).unwrap();
            file.write_all(&raw_len.to_le_bytes()).unwrap();
            file.write_all(&(compressed.len() as u64).to_le_bytes())
                .unwrap();
            file.write_all(&compressed).unwrap();
            directory.push((name, offset, raw_len, compressed.len() as u64));
        }
        file.write_all(&[0]).unwrap();
        let offset = file.stream_position().unwrap();
        file.write_all(&(directory.len() as u32).to_le_bytes())
            .unwrap();
        for (name, position, raw, compressed) in directory {
            file.write_all(&[name.len() as u8]).unwrap();
            file.write_all(name.as_bytes()).unwrap();
            for value in [position, raw, compressed] {
                file.write_all(&value.to_le_bytes()).unwrap();
            }
        }
        file.write_all(&offset.to_le_bytes()).unwrap();
        file.write_all(b"AIED").unwrap();
        drop(file);
        drop(reader);
        let local = load_local_archive("legacy".into(), legacy.clone(), None, false, 0).unwrap();
        let mut collection = assemble_collection(vec![local], true, 1).unwrap();
        std::fs::copy(&legacy, &copy).unwrap();
        collection.archives[0].path = copy.clone();
        collection.archives[0].identity.len = 0;
        let io = validate_archive_identity(&collection.archives[0], false).unwrap();
        assert_eq!(
            io.identity_content_bytes_read,
            std::fs::metadata(&copy).unwrap().len()
        );
        let mut opened = open_source(&collection, 0, None).unwrap();
        assert!(!opened.reader().read("c0").unwrap().is_empty());
        drop(opened);
        collection.archives[0].identity.native_digest = "0".repeat(64);
        assert!(open_source(&collection, 0, None).is_err());
        assert!(validate_archive_identity(&collection.archives[0], false).is_err());
        for path in [rooted, legacy, copy] {
            std::fs::remove_file(path).unwrap();
        }
    }

    #[test]
    fn historical_inode_and_path_collisions_are_not_duplicate_evidence() {
        let base_path = test_path("same-historical-inode-base", "aicollection");
        let child_path = test_path("same-historical-inode-child", "aicollection");
        let base = fixture_collection();
        let mut child = base.clone();
        child.archives[0].id = "B".into();
        child.archives[0].identity.native_digest = "4".repeat(64);
        child.archives[0].identity.encoded_sections_digest = Some("5".repeat(64));
        validate_extension(&base, &child).unwrap();
        write_collection(&base_path, &base).unwrap();
        child.base = Some(BaseCollection {
            path: base_path.clone(),
            root_digest: native_collection_identity(&base_path).unwrap(),
        });
        write_collection(&child_path, &child).unwrap();
        let chain = open_collection_chain(&child_path).unwrap();
        assert_eq!(chain.collection.archives.len(), 2);
        drop(chain);
        child.archives[0].identity.encoded_sections_digest =
            base.archives[0].identity.encoded_sections_digest.clone();
        assert!(validate_extension(&base, &child).is_err());
        for path in [base_path, child_path] {
            std::fs::remove_file(path).unwrap();
        }
    }

    fn hints(path: &Path, entries: serde_json::Value) {
        std::fs::write(
            path,
            serde_json::to_vec(&json!({"schema_version":1,"locations":entries})).unwrap(),
        )
        .unwrap();
    }

    #[test]
    fn copy_ignores_all_historical_stat_fields_and_preserves_io_and_routes() {
        let source = scientific_archive_fixture();
        let copy = test_path("relocated-copy", "aie");
        let (mut collection, _) = routed_collection_fixture(&source);
        let baseline = validate_archive_identity(&collection.archives[0], false).unwrap();
        std::fs::copy(&source, &copy).unwrap();
        let entry = &mut collection.archives[0];
        entry.path = copy.clone();
        entry.identity.len = 0;
        entry.identity.modified_secs = u64::MAX;
        entry.identity.changed_secs = u64::MAX;
        entry.identity.inode = u64::MAX;
        entry.identity.dev = u64::MAX;
        #[cfg(unix)]
        {
            let file = std::fs::OpenOptions::new().write(true).open(&copy).unwrap();
            file.set_times(
                std::fs::FileTimes::new()
                    .set_modified(std::time::UNIX_EPOCH - std::time::Duration::from_secs(1)),
            )
            .unwrap();
        }
        assert_eq!(validate_archive_identity(entry, false).unwrap(), baseline);
        assert!(
            validate_archive_identity(entry, true)
                .unwrap()
                .identity_content_bytes_read
                > 0
        );
        let binding = entry.shape_routes.clone().unwrap();
        let mut opened = open_source(&collection, 0, Some(&binding)).unwrap();
        assert!(!opened.reader().read("shapes").unwrap().is_empty());
        // The relocated file is content-accepted, but a different root is never accepted.
        collection.archives[0].identity.native_digest = "0".repeat(64);
        assert!(open_source(&collection, 0, None).is_err());
        for path in [source, copy] {
            std::fs::remove_file(path).unwrap();
        }
    }

    #[test]
    fn layered_collection_moves_without_rewriting_any_collection_bytes() {
        let source = scientific_archive_fixture();
        let relocated = test_path("layered-archive-copy", "aie");
        let base_path = test_path("relocation-base", "aicollection");
        let moved_base = test_path("relocation-moved-base", "aicollection");
        let child_path = test_path("relocation-child", "aicollection");
        let map = test_path("relocation-locations", "json");
        let (collection, _) = routed_collection_fixture(&source);
        write_collection(&base_path, &collection).unwrap();
        let base_root = native_collection_identity(&base_path).unwrap();
        let mut child = collection.clone();
        child.archives.clear();
        child.junctions.clear();
        child.junction_count = 0;
        child.route_count = 0;
        child.posting_count = 0;
        child.shape_route_blocks.clear();
        child.encoded_shape_route_blocks.clear();
        child.base = Some(BaseCollection {
            path: base_path.clone(),
            root_digest: base_root.clone(),
        });
        write_collection(&child_path, &child).unwrap();
        let child_bytes = std::fs::read(&child_path).unwrap();
        let base_bytes = std::fs::read(&base_path).unwrap();
        let original = open_collection_chain(&child_path).unwrap();
        let baseline = validate_collection_sources(&original.collection, false).unwrap();
        drop(original);
        std::fs::copy(&source, &relocated).unwrap();
        std::fs::rename(&base_path, &moved_base).unwrap();
        std::fs::remove_file(&source).unwrap();
        hints(
            &map,
            json!([
                {"identity":format!("{}:{}",ROOTED_DIRECTORY_SCHEME,collection.archives[0].identity.native_digest),"path":relocated.file_name().unwrap().to_str().unwrap()},
                {"identity":format!("{}:{}",locations::COLLECTION_SCHEME,base_root),"path":moved_base.file_name().unwrap().to_str().unwrap()}
            ]),
        );
        let chain = open_collection_chain_with_locations(&child_path, Some(&map)).unwrap();
        assert_eq!(chain.layers.len(), 2);
        assert_eq!(
            validate_collection_sources(&chain.collection, false).unwrap(),
            baseline
        );
        for layer in &chain.layers {
            let audit = audit_collection_layer(layer).unwrap();
            for routes in audit.shape_routes {
                let a = &layer.manifest.archives[routes.archive_ordinal];
                verify_inspected_archive_routes(
                    a,
                    a.shape_routes.as_ref().unwrap(),
                    &routes.blocks,
                )
                .unwrap();
            }
        }
        assert_eq!(std::fs::read(&child_path).unwrap(), child_bytes);
        assert_eq!(std::fs::read(&moved_base).unwrap(), base_bytes);
        // A location hint cannot replace the committed parent root.
        hints(
            &map,
            json!([{"identity":format!("{}:{}",locations::COLLECTION_SCHEME,base_root),"path":child_path}]),
        );
        assert!(open_collection_chain_with_locations(&child_path, Some(&map)).is_err());
        drop(chain);
        for path in [relocated, moved_base, child_path, map] {
            std::fs::remove_file(path).unwrap();
        }
    }

    #[test]
    fn consumed_corruption_is_not_hidden_by_a_matching_stored_commitment() {
        use std::io::{Seek, SeekFrom};
        let source = scientific_archive_fixture();
        let (collection, _) = routed_collection_fixture(&source);
        let reader = evidence_io::format::SectionReader::open(&source).unwrap();
        let section = reader.section_metadata().find(|s| s.name == "c0").unwrap();
        let offset = section.offset + 1 + section.name.len() as u64 + 16;
        drop(reader);
        let mut file = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .open(&source)
            .unwrap();
        file.seek(SeekFrom::Start(offset)).unwrap();
        let mut byte = [0];
        file.read_exact(&mut byte).unwrap();
        byte[0] ^= 1;
        file.seek(SeekFrom::Start(offset)).unwrap();
        file.write_all(&byte).unwrap();
        drop(file);
        assert_eq!(
            validate_archive_identity(&collection.archives[0], false)
                .unwrap()
                .identity_content_bytes_read,
            0
        );
        assert!(validate_archive_identity(&collection.archives[0], true).is_err());
        let mut archive = open_source(&collection, 0, None).unwrap();
        assert!(archive.reader().read("c0").is_err());
        drop(archive);
        std::fs::remove_file(source).unwrap();
    }

    #[test]
    fn locations_are_strict_bounded_and_relative_to_the_manifest() {
        let map = test_path("strict-locations", "json");
        let id = format!("{}:{}", ROOTED_DIRECTORY_SCHEME, "a".repeat(64));
        hints(&map, json!([{"identity":id,"path":"nested/moved.aie"}]));
        let locations = Locations::load(Some(&map)).unwrap();
        assert_eq!(
            locations.resolve(
                ROOTED_DIRECTORY_SCHEME,
                &"a".repeat(64),
                Path::new("old"),
                Path::new("elsewhere")
            ),
            map.parent().unwrap().join("nested/moved.aie")
        );
        for invalid in [
            json!({"schema_version":2,"locations":[]}),
            json!({"schema_version":1,"locations":[],"unexpected":true}),
            json!({"schema_version":1,"locations":[{"identity":id,"path":"a"},{"identity":id,"path":"b"}]}),
            json!({"schema_version":1,"locations":[{"identity":id,"path":""}]}),
            json!({"schema_version":1,"locations":[{"identity":"unknown:abc","path":"a"}]}),
            json!({"schema_version":1,"locations":[{"identity":format!("{ROOTED_DIRECTORY_SCHEME}:abc"),"path":"a"}]}),
        ] {
            std::fs::write(&map, invalid.to_string()).unwrap();
            assert!(Locations::load(Some(&map)).is_err());
        }
        std::fs::write(&map, vec![b' '; 4 * 1024 * 1024 + 1]).unwrap();
        assert!(Locations::load(Some(&map)).is_err());
        std::fs::remove_file(map).unwrap();
    }
}
