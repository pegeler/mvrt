Version 0.1.0
==============

- Added `mvrt2()` which uses `C++` code to take the place of `mvrtR()`.
- `mvrtR()` now has `max_iterations` parameter to prevent infinite loop
  in cases when very small `max_norm` is specified.
- `mvrtR()` parameters `max.norm` and `type.norm` have been replaced with
  `max_norm` and `type_norm` respectively to align with the transition to
  `mvrt2()`.
- Several addtions to documentation.
