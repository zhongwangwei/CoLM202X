# Material host cleanup phase 5: remove optional scaffolding

Baseline is the verified phase-4 working tree, frozen before edits at
/tmp/colm-host-phase5-_s7lztm7 with full SHA256 manifest. Preserve earlier work.

1. Hoist identical post-hydrology irrigation sum from both MAIN branches to their
   common join; do not move it before a water routine that may modify its inputs.
2. Remove the duplicate thermal diagnostic reset from LandWork preparation.
   THERMAL retains output initialization; standalone PhaseChange retains its own
   reset. Work preparation still owns water snapshots and water-path reset.
3. Factor common parameter/call tails out of MaterialTransport #if/#else pairs,
   preserving exact effective Fortran tokens in MT-on/off configurations.
4. Replace the single history-count forwarding wrapper with its direct material
   call at the generated field entry; update generator and keep all I/O callbacks,
   order, errors, field names and compatibility guards unchanged.

Writer/reviewer separation: parent owns MAIN, LandWork, history generator/wrapper
and integration; independent interface executor owns only signature scaffolding;
independent test owner locks current behavior and compares frozen/current sources.
One final verification batch, with no per-item long simulations. Run targeted raw
bit comparisons for both water branches, thermal preparation/initializer lifecycle
and the generated restart callback; existing phase, snow, thermal and real MPI/I/O
checks; normalized CPP equality for interface-only changes; build changed programs
under supported MT-on/off modes. Use run-7920 for full builds. Since executed
physics and sampling order are unchanged, focused equivalence and compilation
are the primary proof; broaden to long integration only if an unresolved concern
or failure requires it.
