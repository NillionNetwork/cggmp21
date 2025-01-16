Ensures reliability of broadcast channel by adding one extra communication round

CGGMP21 protocol requires a message in the first a round to be sent over reliable
broadcast channel. We ensure reliability of the broadcast channel by introducing extra
communication a round (at cost of additional latency). You may disable it, for instance,
if your transport layer is reliable by construction (e.g. you use blockchain for
communications).

Default: `true`.
