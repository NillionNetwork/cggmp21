/// Refer to `examples/pregenerate.rs` to see how data is pregenerated
#[macro_export]
macro_rules! load_pregenerated_data {
    ($($name:ident: $type:ty),+$(,)?) => {$(
        pub fn $name() -> $type {
            const JSON: &str = include_str!(concat!("../test-data/", stringify!($name), ".json"));
            serde_json::from_str(JSON).unwrap()
        }
    )+};
}
