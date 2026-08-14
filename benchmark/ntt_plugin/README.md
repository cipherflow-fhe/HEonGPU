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

The benchmark executable prints GPUNTT and optimized routes in one run. For
Retina or other HEonGPU application runs, use these runtime routes:

```bash
# Original HEonGPU path: GPUNTT, no PhantomNTT-dependent optimizations.
HEONGPU_USE_PHANTOM_NTT=0 <command>

# Optimized path: PhantomNTT plus KeySwitch_P1, KeySwitch_Part2, and BSGS fusion.
HEONGPU_USE_PHANTOM_NTT=1 <command>
```

With `HEONGPU_USE_PHANTOM_NTT=1`, all KeySwitch/BSGS optimizations are enabled
by default. Use these only for ablation or staged comparison:

```bash
HEONGPU_USE_MOD_KSWITCH=0   # disable KeySwitch_P1
HEONGPU_USE_KSWITCH_P2=0    # disable KeySwitch_Part2
HEONGPU_USE_BSGS_FUSION=0   # disable BSGS fusion
```

The same runtime rule applies to Retina/inference runs: `HEONGPU_USE_PHANTOM_NTT=0`
is the GPUNTT baseline, and `HEONGPU_USE_PHANTOM_NTT=1` is the full optimized
Phantom route unless one of the ablation variables above is explicitly set to
`0`.
