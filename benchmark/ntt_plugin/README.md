# NTT Plugin Benchmarks

Build from `backends/HEonGPU`:

```bash
cmake --build build --target ntt_plugin_benchmark ntt_plugin_correctness --parallel 16
```

Run correctness:

```bash
./build/bin/benchmark/ntt_plugin_correctness
```

Run benchmarks:

```bash
./build/bin/benchmark/ntt_plugin_benchmark
```

Use `HEONGPU_USE_PHANTOM_NTT=0` to force the GPUNTT path.
Use `HEONGPU_USE_PHANTOM_NTT=1` to enable PhantomNTT; KeySwitch_P1 is enabled
by default in this mode.

For comparison only, `HEONGPU_USE_PHANTOM_NTT=1 HEONGPU_USE_MOD_KSWITCH=0`
keeps PhantomNTT on and disables KeySwitch_P1.
