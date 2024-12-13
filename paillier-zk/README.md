![License](https://img.shields.io/crates/l/paillier-zk.svg)
[![Docs](https://docs.rs/paillier-zk/badge.svg)](https://docs.rs/paillier-zk)
[![Crates io](https://img.shields.io/crates/v/paillier-zk.svg)](https://crates.io/crates/paillier-zk)
[![Discord](https://img.shields.io/discord/905194001349627914?logo=discord&logoColor=ffffff&label=Discord)][in Discord]

[in Discord]: https://discordapp.com/channels/905194001349627914/1285268686147424388

# paillier-zk

This crate provides ZK-proofs for some properties about paillier encryption.
See the module docs for the properties and examples of usage.

This library is built on top of [fast-paillier](https://lib.rs/fast-paillier) crate.
This crate and the underlying big integer implementation are reexported for the
consumer to be able to use them, instead of trying to match a version.
