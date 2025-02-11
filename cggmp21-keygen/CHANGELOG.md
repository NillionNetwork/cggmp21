# Changelog

## v0.5.0
* Update `hd-wallet` dep to v0.6 [#120]

[#120]: https://github.com/LFDT-Lockness/cggmp21/pull/120

## v0.4.0
* BREAKING: use `hd-wallet` crate for HD support instead of `slip-10` [#115]
* BREAKING: rename `hd-wallets` feature into `hd-wallet` [#115]
* Update `key-share` to v0.5

[#115]: https://github.com/LFDT-Lockness/cggmp21/pull/115

## v0.3.2
* Update links in the documentation and crate settings after moving the repo [#113]

[#113]: https://github.com/LFDT-Lockness/cggmp21/pull/113

## v0.3.1
* Take advantage of `#[udigest(as = ...)]` attribute [#106]

[#106]: https://github.com/LFDT-Lockness/cggmp21/pull/106

## v0.3.0
* security fix: derive challenges for zero-knowledge proof unambiguously
* Update `udigest` to v0.2
* Update `generic-ec` to v0.4
* Update `slip-10` to v0.4

## v0.2.0
* Make library `#![no_std]`-compatible and WASM-friendly [#100]
* Provide sync API to carry out DKG protocol [#100]
* Update `round-based` dep to `v0.3`
* Update `generic-ec` and `slip-10` deps to latest version [#101]

[#100]: https://github.com/LFDT-Lockness/cggmp21/pull/100
[#101]: https://github.com/LFDT-Lockness/cggmp21/pull/101

## v0.1.0

Initial release
