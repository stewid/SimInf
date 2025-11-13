# Class `"SimInf_individual_events"`

Class `"SimInf_individual_events"`

## Slots

- `id`:

  an integer or character identifier of the individual.

- `event`:

  four event types are supported: *exit*, *enter*, *internal transfer*,
  and *external transfer*. When assigning the events, they can either be
  coded as a numerical value or a character string: *exit;* `0` or
  `'exit'`, *enter;* `1` or `'enter'`, *internal transfer;* `2` or
  `'intTrans'`, and *external transfer;* `3` or `'extTrans'`.

- `time`:

  an integer, character, or date (of class `Date`) for when the event
  occured. If it's a character it must be able to coerce to `Date`.

- `node`:

  an integer or character identifier of the source node.

- `dest`:

  an integer or character identifier of the destination node.
