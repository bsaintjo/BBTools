#!/bin/bash

# simd_incubator_tools.sh
#
# Identifies BBTools bash scripts whose Java code paths will instantiate
# objects from jdk.incubator.vector, and explains exactly how each one
# reaches that code — including which scripts do so without simd=t.
#
# ─── BACKGROUND: HOW $SIMD AND Shared.SIMD RELATE ───────────────────────────
#
# There are TWO separate "simd" concepts in BBTools:
#
#   $SIMD  (shell variable)
#       = "--add-modules jdk.incubator.vector"
#       Set by javasetup.sh / calcmem.sh on the JVM command line.
#       Controls whether the jdk.incubator.vector module is accessible at all.
#       Auto-detected (no user flag needed) when AVX2/AVX-512/SVE + Java 17+.
#       Can be forced on/off with simd=t / simd=f on the shell command line.
#
#   Shared.SIMD  (Java boolean)
#       Controls whether guarded code paths (if(Shared.SIMD){...}) execute.
#       Set by: Vector.simd256  (auto, at class-load time)
#               OR by the Java "simd=t" argument parsed from args[]
#               OR hardcoded in a few entry-point main() methods.
#
# ─── THE STATIC-INIT CHAIN (fires without simd=t) ────────────────────────────
#
# When --add-modules jdk.incubator.vector is on the JVM command line, the
# following chain fires at startup for every tool that references shared.Shared:
#
#   shared.Shared (loaded by virtually every BBTools class)
#     └─ Shared.java:132   SIMD = (Vector.simd256)
#          └─ loads simd.Vector
#               ├─ Vector.java:1222  vectorLoaded = vectorLoaded()
#               │    └─ Class.forName("jdk.incubator.vector.ByteVector")
#               │         → ByteVector class loaded      ← INCUBATOR OBJECT
#               └─ Vector.java:1223  maxSimdWidth = maxSimdWidth()
#                    └─ SIMD.maxVectorLength()
#                         → loads simd.SIMD
#                              ├─ SIMD.java:24  ByteVector.SPECIES_PREFERRED  ← INCUBATOR OBJECT
#                              ├─ SIMD.java:28  FloatVector.SPECIES_256       ← INCUBATOR OBJECT
#                              ├─ SIMD.java:33  IntVector.SPECIES_256         ← INCUBATOR OBJECT
#                              ├─ SIMD.java:38  ShortVector.SPECIES_256       ← INCUBATOR OBJECT
#                              ├─ SIMD.java:43  DoubleVector.SPECIES_256      ← INCUBATOR OBJECT
#                              └─ SIMD.java:48  LongVector.SPECIES_256        ← INCUBATOR OBJECT
#
# Result: Shared.SIMD = (maxSimdWidth >= 256), i.e. true on AVX2+ hardware.
# No user-supplied simd=t is needed for this chain to fire.
#
# ─── SCRIPTS THAT ALSO HARDCODE Shared.SIMD=true IN JAVA ────────────────────
#
# aligner/AlignRandom.java:38 sets Shared.SIMD=true unconditionally in main(),
# before any argument parsing.  This is the only entry-point class that does so.
# It is launched by alignrandom.sh.
#
# ─── SCRIPTS THAT LOAD ADDITIONAL INCUBATOR CLASSES ─────────────────────────
#
# Beyond simd.SIMD, the following classes also have jdk.incubator.vector
# VectorSpecies objects as static fields (instantiated at class-load time):
#
#   simd/SIMDByte256.java:25-31
#       ByteVector.SPECIES_256, ByteVector.SPECIES_64, LongVector.SPECIES_64,
#       plus a static{} block building VectorShuffle objects.
#       Loaded when Shared.SIMD=true and simd.Vector dispatches to SIMDByte256.
#       → Affects all read-processing tools (reformat, bbduk, bbmap, etc.)
#
#   simd/SIMDAlign.java:25-55
#       FloatVector/ByteVector/IntVector/ShortVector/DoubleVector/LongVector SPECIES_256.
#       Loaded when idaligner.* classes call SIMDAlign methods.
#       → Affects all idaligner.* tools.
#
#   simd/SIMDAlignByte.java:28-46
#       Same species set as SIMDAlign.
#       Loaded when ifa.IndelFreeAligner* or idaligner.IDAlignerStatics runs.
#       → Affects indelfree.sh and all idaligner.* tools.
#
#   idaligner/DiagonalAligner.java:26-27
#       ByteVector.SPECIES_128, IntVector.SPECIES_128 — static fields.
#       Loaded unconditionally when idaligner.Test or TestAlignerSuite runs.
#       → Affects testaligners.sh, testaligners2.sh.
#
#   rand/FastRandomSIMD.java:23
#       LongVector.SPECIES_256 — static field.
#       Loaded unconditionally when FastRandomSIMD is instantiated.
#       No dedicated bash script; used internally.
#
#   fun/ByteVectorBench.java:18-20
#       ByteVector.SPECIES_64/128/256 — static fields.
#       No dedicated bash script; run directly.

echo "========================================================================"
echo " BBTools scripts that instantiate jdk.incubator.vector objects"
echo "========================================================================"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
echo "── alignrandom.sh ───────────────────────────────────────────────────────"
echo "   Java class:  aligner.AlignRandom"
echo "   Why special: AlignRandom.main() sets Shared.SIMD=true UNCONDITIONALLY"
echo "                (AlignRandom.java:38) before any argument parsing."
echo "                It then uses GlocalPlusAligner5, which calls"
echo "                simd.SIMDAlign.alignBandVector() — loading SIMDAlign's"
echo "                six VectorSpecies static fields."
echo "   Incubator objects created regardless of CPU, Java version, or flags."
echo "   Requires: --add-modules jdk.incubator.vector  (set by \$SIMD auto-detect"
echo "             or explicitly; without it the module is inaccessible and"
echo "             Class.forName fails silently in vectorLoaded())."
echo ""

# ─────────────────────────────────────────────────────────────────────────────
echo "── indelfree.sh ─────────────────────────────────────────────────────────"
echo "   Java class:  ifa.IndelFreeAligner4"
echo "   Path to incubator:"
echo "     1. Shared.SIMD = Vector.simd256  (auto-detected via static-init chain)"
echo "     2. IndelFreeAligner4.alignAllPositions() checks Shared.SIMD"
echo "        → calls simd.SIMDAlignByte.alignDiagonal()  (IndelFreeAligner4.java:452)"
echo "        → SIMDAlignByte loaded: ByteVector/FloatVector/IntVector/ShortVector/"
echo "          LongVector SPECIES_256 static fields instantiated"
echo "   Note: 'simd=t' in the usage text is the documented default behavior;"
echo "         the Java class does NOT hardcode simd=t — it relies on auto-detect."
echo ""

# ─────────────────────────────────────────────────────────────────────────────
echo "── testaligners.sh / testaligners2.sh ───────────────────────────────────"
echo "   Java classes: idaligner.Test / idaligner.TestAlignerSuite"
echo "   Path to incubator:"
echo "     1. Shared.SIMD = Vector.simd256  (static-init chain)"
echo "     2. Test/TestAlignerSuite instantiates DiagonalAligner"
echo "        → DiagonalAligner loaded: ByteVector.SPECIES_128 and"
echo "          IntVector.SPECIES_128 static fields instantiated (DiagonalAligner.java:26-27)"
echo "          — these fire at class-load time, no Shared.SIMD guard."
echo "     3. Also loads SIMDAlignByte and SIMDAlign via other IDAligner impls."
echo ""

# ─────────────────────────────────────────────────────────────────────────────
echo "── idaligner scripts (bandedaligner, bandedplusaligner, crosscutaligner,"
echo "   driftingaligner, driftingplusaligner, glocalaligner, quantumaligner,"
echo "   scrabblealigner, wavefrontaligner, wavefrontalignerviz, wobblealigner,"
echo "   wobbleplusaligner, quabblealigner, xdrophaligner) ─────────────────────"
echo "   Java classes: idaligner.BandedAligner, BandedPlusAligner2,"
echo "                 CrossCutAligner, DriftingAligner, DriftingPlusAligner2,"
echo "                 GlocalAligner, QuantumAligner, ScrabbleAligner2,"
echo "                 WaveFrontAligner, WaveFrontAlignerViz, WobbleAligner,"
echo "                 WobblePlusAligner3, QuabbleAligner, XDropHAligner"
echo "   Path to incubator:"
echo "     1. Shared.SIMD = Vector.simd256  (static-init chain)"
echo "     2. Each aligner's inner loop checks Shared.SIMD"
echo "        → calls simd.SIMDAlign.alignBandVector() (or variant)"
echo "        → SIMDAlign loaded: FloatVector/ByteVector/IntVector/ShortVector/"
echo "          DoubleVector/LongVector SPECIES_256 static fields instantiated"
echo "     3. Also calls simd.SIMDAlignByte via idaligner.IDAlignerStatics"
echo "        → SIMDAlignByte loaded with same species set"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
echo "── Read-processing tools (reformat, bbduk, bbmap, bbmerge, bbnorm, seal,"
echo "   clumpify, filescan, fastqscan, stream, pileup, callvariants, sketch,"
echo "   rqcfilter3, tadpole, bbcms, kmercountexact, stats, callgenes,"
echo "   consensus, icecreamfinder, train, seqtovec, scalars, cloudplot, ...)"
echo "   ────────────────────────────────────────────────────────────────────────"
echo "   Path to incubator:"
echo "     1. Shared.SIMD = Vector.simd256  (static-init chain)"
echo "     2. Read processing calls simd.Vector methods (countMatches, sum,"
echo "        reverseComplementInPlace, findSymbols, etc.)"
echo "        → simd.Vector dispatches to SIMDByte256 when Shared.SIMD=true"
echo "        → SIMDByte256 loaded: ByteVector.SPECIES_256, ByteVector.SPECIES_64,"
echo "          LongVector.SPECIES_64 static fields + VectorShuffle static block"
echo "   These tools do NOT need simd=t; auto-detection on AVX2 + Java 17+ suffices."
echo ""

# ─────────────────────────────────────────────────────────────────────────────
echo "========================================================================"
echo " Summary: which scripts need simd=t vs. which don't"
echo "========================================================================"
echo ""
echo "  alignrandom.sh      — Shared.SIMD=true hardcoded in Java main()."
echo "                        Incubator objects created unconditionally"
echo "                        (as long as --add-modules is on the JVM line)."
echo ""
echo "  All other scripts   — Rely on auto-detection:"
echo "                        $SIMD set by javasetup.sh when AVX2 + Java 17+."
echo "                        Shared.SIMD = Vector.simd256 via static-init chain."
echo "                        No simd=t needed from the user."
echo ""
echo "  simd=t (user flag)  — Forces \$SIMD and Shared.SIMD=true even on"
echo "                        non-AVX2 hardware or Java < 17."
echo "                        No bash script hardcodes simd=t in its launch cmd."
echo ""
echo "  simd=f (user flag)  — Disables \$SIMD and prevents the module from"
echo "                        loading; also sets Shared.SIMD=false."
echo "                        Prevents ALL incubator object instantiation."
echo ""
