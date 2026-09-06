//! Allocator selection belongs to the executable, never to the reusable libraries.
//! Do not enable mimalloc's `override`: native dependencies retain their own allocator.
//! Cargo selects the measured v2 backend; the default build remains the system allocator.
#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

pub(crate) const fn name() -> &'static str {
    if cfg!(feature = "mimalloc") {
        "mimalloc"
    } else {
        "system"
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn aligned_growth_and_cross_thread_free() {
        #[repr(align(4096))]
        struct Aligned([u8; 8192]);
        let aligned = Box::new(Aligned([17; 8192]));
        assert_eq!(aligned.0.as_ptr() as usize % 4096, 0);
        let mut data = vec![3u64; 17];
        data.resize(100_000, 42);
        let data = std::thread::spawn(move || {
            assert_eq!(aligned.0[8191], 17);
            drop(aligned);
            assert!(data[..17].iter().all(|&v| v == 3));
            assert!(data[17..].iter().all(|&v| v == 42));
            data.truncate(50);
            data.shrink_to_fit();
            data
        })
        .join()
        .unwrap();
        assert_eq!(data.len(), 50);
        assert_eq!(data[49], 42);
    }
}
