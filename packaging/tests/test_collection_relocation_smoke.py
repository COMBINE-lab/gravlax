import importlib.util
import os
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location(
    'collection_relocation_smoke',
    Path(__file__).resolve().parents[1] / 'collection_relocation_smoke.py',
)
SMOKE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SMOKE)


class RelocationComparisonTests(unittest.TestCase):
    def test_macos_canonical_and_manifest_paths_match(self):
        aliases = ('/private/var/tmp/original/bundle', '/var/tmp/original/bundle',
                   '/private/var/tmp/copy/bundle', '/var/tmp/copy/bundle')
        before = {'path': aliases[0] + '/a.aie'}
        after = {'path': aliases[3] + '/a.aie'}
        self.assertEqual(SMOKE.normalized_result(before, aliases),
                         SMOKE.normalized_result(after, aliases))

    def test_windows_verbatim_drive_paths_match(self):
        aliases = SMOKE.bundle_path_aliases(Path(r'C:\Temp\original'), Path(r'C:\Temp\copy'))
        self.assertEqual(SMOKE.normalized_result(r'\\?\C:\Temp\original\a.aie', aliases),
                         SMOKE.normalized_result(r'C:\Temp\copy\a.aie', aliases))

    def test_windows_verbatim_unc_paths_match(self):
        aliases = SMOKE.bundle_path_aliases(Path(r'\\server\share\original'),
                                          Path(r'\\server\share\copy'))
        self.assertEqual(SMOKE.normalized_result(r'\\?\UNC\server\share\original\a.aie', aliases),
                         SMOKE.normalized_result(r'\\server\share\copy\a.aie', aliases))

    @unittest.skipIf(os.name == 'nt', 'symlink creation can require Windows privileges')
    def test_canonical_alias_is_captured_before_original_disappears(self):
        with tempfile.TemporaryDirectory() as temporary:
            parent = Path(temporary)
            actual = parent / 'actual'
            actual.mkdir()
            alias = parent / 'alias'
            alias.symlink_to(actual, target_is_directory=True)
            original = alias / 'bundle'
            original.mkdir()
            paths = SMOKE.bundle_path_aliases(original)
            original.rename(actual / 'unavailable')
            alias.unlink()
            self.assertEqual(SMOKE.normalized_result(str(actual.resolve() / 'bundle/a.aie'), paths),
                             '<bundle>/a.aie')

    def test_only_complete_bundle_prefixes_are_normalized(self):
        aliases = ('/tmp/bundle',)
        for value in ('/tmp/bundle-other/a.aie', '/elsewhere/tmp/bundle/a.aie',
                      'sample /tmp/bundle/a.aie', 'cdb2cc2b3a79b388'):
            self.assertEqual(SMOKE.normalized_result(value, aliases), value)
        self.assertEqual(SMOKE.normalized_result('/tmp/bundle', aliases), '<bundle>')

    def test_content_io_and_filename_differences_remain_visible(self):
        before = {'path': '/tmp/original/a.aie', 'root_digest': 'abc',
                  'io': {'source_identity_total_bytes_read': 3993},
                  'rows': [[1, 2]], 'elapsed_seconds': 1.0}
        aliases = ('/tmp/original', '/tmp/copy')
        after = {**before, 'path': '/tmp/copy/a.aie', 'elapsed_seconds': 2.0,
                 'locations_manifest': {'path': '/tmp/copy/locations.json'}}
        expected = SMOKE.normalized_result(before, aliases)
        self.assertEqual(expected, SMOKE.normalized_result(after, aliases))
        for changed in ({'path': '/tmp/copy/b.aie'}, {'root_digest': 'def'},
                        {'io': {'source_identity_total_bytes_read': 3994}},
                        {'rows': [[1, 3]]}):
            with self.subTest(changed=changed):
                self.assertNotEqual(expected, SMOKE.normalized_result({**after, **changed}, aliases))


if __name__ == '__main__':
    unittest.main()
