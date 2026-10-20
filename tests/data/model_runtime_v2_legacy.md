The companion `model_runtime_v2_legacy.bin` was generated using the runtime,
implicit rod, IB adapter, and force-coupled source at commit aa91785, before the
ModelRuntime/IBCoupling split. It uses the existing native version-2 checkpoint
layout on little-endian x86-64 with double-precision real_t.

The fixture contains one BendingRodImplicit (persistent type 14), default
configuration except links=2, advanced once at dt=0.01 with zero applied load.
It includes generalized-alpha acceleration history and IB force storage.
The regression restores it and compares the state and next step against direct
construction. This format was already platform-native; the fixture does not
introduce a portable checkpoint format.
