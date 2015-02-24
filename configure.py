#!/usr/bin/python

from __future__ import print_function
import os
import sys
import argparse
import ninja


def get_program_info(program, arguments, use_stdout=True):
    from subprocess import PIPE, Popen

    if not isinstance(arguments, list):
        arguments = [str(arguments)]
    process = Popen([program] + arguments, stdout=PIPE, stderr=PIPE, bufsize=1)
    outdata, errdata = process.communicate()
    if use_stdout:
        return outdata
    else:
        return errdata


def detect_intel_compiler(program):
    banner = get_program_info(program, "-v", use_stdout=False)
    if banner:
        version = banner.splitlines()[0]
        import re
        return bool(re.match("(icc|ifort) version (\d+(:?\.\d+)+)", version))
    return False


class Configuration:
    _fflags_intel_map = {
        "penryn": ["-m64", "-xSSE4.1"],
        "nehalem": ["-m64", "-xSSE4.2"],
        "sandybridge": ["-m64", "-xAVX"],
        "ivybridge": ["-m64", "-xCORE-AVX-I"],
        "haswell": ["-m64", "-xCORE-AVX2"],
        "mic": ["-mmic"],
        "piledriver": ["-m64", "-mavx", "-fma"]
    }
    _fflags_gnu_map = {
        "penryn": ["-m64", "-march=core2", "-msse4.1"],
        "nehalem": ["-m64", "-march=corei7"],
        "sandybridge": ["-m64", "-march=corei7-avx"],
        "ivybridge": ["-m64", "-march=core-avx-i"],
        "haswell": ["-m64", "-march=core-avx2"],
        "piledriver": ["-m64", "-march=bdver2"]
    }
    _cflags_intel_map = {
        "penryn": ["-m64", "-xSSE4.1", "-DERD_PNR"],
        "nehalem": ["-m64", "-xSSE4.2", "-DERD_NHM"],
        "sandybridge": ["-m64", "-xAVX", "-DERD_SNB"],
        "ivybridge": ["-m64", "-xCORE-AVX-I", "-DERD_IVB"],
        "haswell": ["-m64", "-xCORE-AVX2", "-DERD_HSW"],
        "mic": ["-mmic", "-no-opt-prefetch", "-DEFD_MIC"],
        "piledriver": ["-m64", "-mavx", "-fma", "-DERD_PLD"]
    }
    _cflags_gnu_map = {
        "penryn": ["-m64", "-march=core2", "-msse4.1", "-DERD_PNR"],
        "nehalem": ["-m64", "-march=corei7", "-DERD_NHM"],
        "sandybridge": ["-m64", "-march=corei7-avx", "-DERD_SNB"],
        "ivybridge": ["-m64", "-march=core-avx-i", "-DERD_IVB"],
        "haswell": ["-m64", "-march=core-avx2", "-DERD_HSW"],
        "piledriver": ["-m64", "-march=bdver2", "-DERD_PLD"]
    }
    _extra_fflags = ["-O3", "-g", "-reentrancy", "threaded", "-recursive"]
    _extra_cflags = ["-O3", "-g", "-std=gnu99", "-D__ALIGNLEN__=64", "-Wall", "-Wextra", "-Werror", "-Wno-unused-variable", "-openmp"]
    _native_cflags = ["-D__ERD_PROFILE__"]
    _offload_cflags = ["-offload-option,mic,compiler,\"-z defs -no-opt-prefetch\""]
    _extra_ldflags = ["-static-intel", "-no-intel-extensions", "-lifcore", "-lrt", "-openmp"]

    def __init__(self, options, ninja_build_file=os.path.join(os.path.dirname(os.path.abspath(__file__)), "build.ninja")):
        self.output = open(ninja_build_file, "w")
        self.writer = ninja.Writer(self.output)
        self.source_dir = None
        self.build_dir = None
        self.artifact_dir = None
        self.object_ext = ".o"

        use_icc = detect_intel_compiler(options.cc)
        use_ifort = detect_intel_compiler(options.fc)

        if use_icc:
            cflags = Configuration._cflags_intel_map[options.arch]
        else:
            cflags = Configuration._cflags_gnu_map[options.arch]
        if options.offload:
            cflags += Configuration._offload_cflags
        else:
            cflags += ["-D__ERD_PROFILE__"]
            if use_icc:
                cflags += ["-offload=none", "-diag-disable", "161,2423"]
        cflags += Configuration._extra_cflags
        fflags = []
        if use_ifort:
            fflags = Configuration._fflags_intel_map[options.arch]
        else:
            fflags = Configuration._fflags_gnu_map[options.arch]
        fflags += Configuration._extra_fflags
        ldflags = Configuration._extra_ldflags

        cflags = " ".join(cflags)
        if options.cflags:
            cflags = cflags + " " + options.cflags
        fflags = " ".join(fflags)
        if options.fflags:
            fflags = fflags + " " + options.fflags
        ldflags = " ".join(ldflags)
        if options.ldflags:
            ldflags = ldflags + " " + options.ldflags

        # Variables
        self.writer.variable("fc", options.fc)
        self.writer.variable("fflags", fflags)
        self.writer.variable("cc", options.cc)
        self.writer.variable("cflags", cflags)
        self.writer.variable("ldflags", ldflags)
        self.writer.variable("ar", options.ar)
        if options.offload:
            self.writer.variable("arflags", "-qoffload-build")

        # Rules
        self.writer.rule("CC", "$cc $cflags $includes -MMD -MT $out -MF $out.d -o $out -c $in",
            description="CC $descpath",
            depfile="$out.d")
        self.writer.rule("FC", "$fc $fflags -o $out -c $in",
            description="FC $descpath")
        self.writer.rule("CCLD", "$cc -o $out $objs $libdirs $libs $ldflags",
            description="CCLD $descpath")
        self.writer.rule("AR", "$ar $arflags rcs $out $in",
            description="AR $descpath")

    def cc(self, source_file, object_file=None, include_dirs=[]):
        if not os.path.isabs(source_file):
            source_file = os.path.join(self.source_dir, source_file)
        if object_file is None:
            object_file = os.path.join(self.build_dir, os.path.relpath(source_file, self.source_dir)) + self.object_ext
        variables = {
            "descpath": os.path.relpath(source_file, self.source_dir)
        }
        if include_dirs:
            variables["includes"] = " ".join(["-I" + i for i in include_dirs])
        self.writer.build(object_file, "CC", source_file, variables=variables)
        return object_file

    def fc(self, source_file, object_file=None):
        if not os.path.isabs(source_file):
            source_file = os.path.join(self.source_dir, source_file)
        if object_file is None:
            object_file = os.path.join(self.build_dir, os.path.relpath(source_file, self.source_dir)) + self.object_ext
        self.writer.build(object_file, "FC", source_file,
            variables={"descpath": os.path.relpath(source_file, self.source_dir)})
        return object_file

    def ccld(self, object_files, executable_file, lib_dirs=[], libs=[], lib_deps=[]):
        if not os.path.isabs(executable_file):
            executable_file = os.path.join(self.artifact_dir, executable_file)
        variables = {
            "descpath": os.path.relpath(executable_file, self.artifact_dir)
        }
        if lib_dirs:
            variables["libdirs"] = " ".join(["-L" + l for l in lib_dirs])
        if libs:
            variables["libs"] = " ".join(["-l" + l for l in libs])
        variables["objs"] = " ".join(object_files)
        self.writer.build(executable_file, "CCLD", object_files + lib_deps, variables=variables)
        return executable_file

    def ar(self, object_files, archive_file):
        if not os.path.isabs(archive_file):
            archive_file = os.path.join(self.artifact_dir, archive_file)
        self.writer.build(archive_file, "AR", object_files,
            variables={"descpath": os.path.relpath(archive_file, self.artifact_dir)})
        return archive_file


parser = argparse.ArgumentParser(
    description="OptErd configuration script")
parser.add_argument("-arch", dest="arch", required=True,
    choices=("penryn", "nehalem", "sandybridge", "ivybridge", "haswell", "mic", "piledriver"),
    help="Target microarchitecture")
parser.add_argument("-enable-offload", dest="offload", action="store_true",
    help="Enable MIC offload")
parser.add_argument("--with-cc", dest="cc", default=os.getenv("CC", "icc"))
parser.add_argument("--with-cflags", dest="cflags", default=os.getenv("CFLAGS"))
parser.add_argument("--with-fc", dest="fc", default=os.getenv("FC", "ifort"))
parser.add_argument("--with-fflags", dest="fflags", default=os.getenv("FFLAGS"))
parser.add_argument("--with-ar", dest="ar", default=os.getenv("AR", "xiar"))
parser.add_argument("--with-ldflags", dest="ldflags", default=os.getenv("LDFLAGS"))


erd_sources = [
    "erd__memory_csgto.c",
    "erd__1111_csgto.c", "erd__2d_coefficients.c", "erd__2d_pq_integrals.c",
    "erd__boys_table.c", "erd__jacobi_table.c", "erd__cartesian_norms.c", "erd__csgto.c",
    "erd__dsqmin_line_segments.c", "erd__e0f0_pcgto_block.c", "erd__hrr_matrix.c",
    "erd__hrr_step.c", "erd__hrr_transform.c", "erd__int2d_to_e000.c", "erd__int2d_to_e0f0.c",
    "erd__move_ry.c", "erd__normalize_cartesian.c",
    "erd__pppp_pcgto_block.c", "erd__rys_1_roots_weights.c", "erd__rys_2_roots_weights.c", "erd__rys_3_roots_weights.c",
    "erd__rys_4_roots_weights.c", "erd__rys_5_roots_weights.c", "erd__rys_roots_weights.c", "erd__rys_x_roots_weights.c",
    "erd__set_abcd.c", "erd__set_ij_kl_pairs.c", "erd__spherical_transform.c", "erd__sppp_pcgto_block.c",
    "erd__sspp_pcgto_block.c", "erd__sssp_pcgto_block.c", "erd__ssss_pcgto_block.c", "erd__xyz_to_ry_abcd.c",
    "erd__xyz_to_ry_matrix.c",
    "erd_profile.c"
]

oed_sources = [
    "oed__cartesian_norms.f", "oed__ctr_2index_block.f", "oed__ctr_2index_reorder.f", "oed__ctr_3index_block.f", "oed__ctr_3index_reorder.f",
    "oed__ctr_pair_new.f", "oed__ctr_pair_update.f", "oed__ctr_rs_expand.f", "oed__ctr_single_new.f", "oed__ctr_single_update.f",
    "oed__dsqmin_line_segments.f", "oed__gener_kin_batch.f", "oed__gener_kin_derv_batch.f", "oed__gener_nai_batch.f", "oed__gener_nai_derv_batch.f",
    "oed__gener_ovl3c_batch.f", "oed__gener_ovl_batch.f", "oed__gener_ovl_derv_batch.f", "oed__gener_xyz_batch.f", "oed__gener_xyz_derv_batch.f",
    "oed__hrr_matrix.f", "oed__hrr_step.f", "oed__hrr_transform.f", "oed__kin_1d_derv_integrals.f", "oed__kin_1d_integrals.f",
    "oed__kin_ab_def_blocks.f", "oed__kin_ab_pcgto_block.f", "oed__kin_csgto.f", "oed__kin_derv_csgto.f", "oed__kin_derv_def_blocks.f",
    "oed__kin_derv_int1d_to_00.f", "oed__kin_derv_int1d_to_a0.f", "oed__kin_derv_int1d_to_ab.f", "oed__kin_derv_pcgto_block.f", "oed__kin_int1d_to_a0.f",
    "oed__kin_int1d_to_ab.f", "oed__kin_prepare_ctr.f", "oed__kin_set_ab.f", "oed__kin_set_derv_ab.f", "oed__kin_set_derv_sequence.f",
    "oed__kin_set_ij_pairs.f", "oed__map_ijkl_to_ikjl.f", "oed__memory_kin_batch.f", "oed__memory_kin_csgto.f", "oed__memory_kin_derv_batch.f",
    "oed__memory_kin_derv_csgto.f", "oed__memory_nai_batch.f", "oed__memory_nai_csgto.f", "oed__memory_nai_derv_batch.f", "oed__memory_nai_derv_csgto.f",
    "oed__memory_ovl3c_batch.f", "oed__memory_ovl3c_csgto.f", "oed__memory_ovl_batch.f", "oed__memory_ovl_csgto.f", "oed__memory_ovl_derv_batch.f",
    "oed__memory_ovl_derv_csgto.f", "oed__move_ry.f", "oed__nai_1d_ab_integrals.f", "oed__nai_1d_cenderv_integrals.f", "oed__nai_1d_coefficients.f",
    "oed__nai_1d_p_integrals.f", "oed__nai_1d_shderv_integrals.f", "oed__nai_csgto.f", "oed__nai_derv_2cen_pcgto_block.f", "oed__nai_derv_3cen_pcgto_block.f",
    "oed__nai_derv_csgto.f", "oed__nai_derv_def_blocks.f", "oed__nai_derv_int1d_to_00.f", "oed__nai_derv_int1d_to_a0.f", "oed__nai_derv_int1d_to_ab.f",
    "oed__nai_e0_def_blocks.f", "oed__nai_e0_pcgto_block.f", "oed__nai_int1d_to_e0.f", "oed__nai_prepare_ctr.F", "oed__nai_set_ab.f",
    "oed__nai_set_derv_ab.f", "oed__nai_set_derv_ijc_triples.f", "oed__nai_set_ijc_triples.f", "oed__normalize_cartesian.f", "oed__ovl_1d_ab_integrals.f",
    "oed__ovl_1d_derv_integrals.f", "oed__ovl_1d_integrals.f", "oed__ovl3c_csgto.f", "oed__ovl3c_f00_def_blocks.f", "oed__ovl3c_f00_pcgto_block.f",
    "oed__ovl3c_prepare_ctr.F", "oed__ovl3c_set_abc.f", "oed__ovl3c_set_ijk.f", "oed__ovl3c_set_ijk_triples.f", "oed__ovl_csgto.f",
    "oed__ovl_derv_csgto.f", "oed__ovl_derv_def_blocks.f", "oed__ovl_derv_int1d_to_00.f", "oed__ovl_derv_int1d_to_a0.f", "oed__ovl_derv_int1d_to_ab.f",
    "oed__ovl_derv_pcgto_block.f", "oed__ovl_e0_def_blocks.f", "oed__ovl_e0_pcgto_block.f", "oed__ovl_int1d_to_e0.f", "oed__ovl_prepare_ctr.F",
    "oed__ovl_set_ab.f", "oed__ovl_set_derv_ab.f", "oed__ovl_set_derv_sequence.f", "oed__ovl_set_ij_pairs.f", "oed__print_2ind_batch.f",
    "oed__print_batch.f", "oed__rys_1_roots_weights.f", "oed__rys_2_roots_weights.f", "oed__rys_3_roots_weights.f", "oed__rys_4_roots_weights.f",
    "oed__rys_5_roots_weights.f", "oed__rys_roots_weights.f", "oed__rys_x_roots_weights.f", "oed__spherical_transform.f", "oed__transpose_batch.f",
    "oed__xyz_1d_ab_integrals.f", "oed__xyz_1d_derv_integrals.f", "oed__xyz_1d_integrals.f", "oed__xyz_1d_mom_ab_integrals.f",
    "oed__xyz_1d_ovl_mom_ab_integrals.f", "oed__xyz_csgto.f", "oed__xyz_derv_csgto.f", "oed__xyz_derv_def_blocks.f", "oed__xyz_derv_int1d_to_00.f",
    "oed__xyz_derv_int1d_to_a0.f", "oed__xyz_derv_int1d_to_ab.f", "oed__xyz_derv_pcgto_block.f", "oed__xyz_e0_def_blocks.f", "oed__xyz_e0_pcgto_block.f",
    "oed__xyz_init_00_ints.f", "oed__xyz_init_ovrlp.f", "oed__xyz_int1d_to_e0.f", "oed__xyz_mom_integrals.f", "oed__xyz_set_ab.f", "oed__xyz_set_derv_ab.f",
    "oed__xyz_set_derv_sequence.f", "oed__xyz_to_ry_abc.f", "oed__xyz_to_ry_ab.f", "oed__xyz_to_ry_matrix.f"
]

cint_sources = ["cint_basisset.c", "erd_integral.c", "oed_integral.c", "cint_offload.c"]


def main():
    options = parser.parse_args()

    config = Configuration(options)


    # print('rule GENERATE_HEADER', file = makefile)
    # print(tab + 'command = cd $WORKDIR && bash $SCRIPT', file = makefile)
    # print(tab + 'description = SH $SCRIPT', file = makefile)

    root_dir = os.path.dirname(os.path.abspath(__file__))

    # Build ERD library
    config.source_dir = os.path.join(root_dir, "erd")
    config.build_dir = os.path.join(root_dir, "build", "erd")
    config.artifact_dir = os.path.join(root_dir, "lib")
    erd_objects = []
    for erd_source in erd_sources:
        erd_objects.append(config.cc(erd_source))
    liberd = config.ar(erd_objects, "liberd.a")

    # Build OED library
    config.source_dir = os.path.join(root_dir, "oed")
    config.build_dir = os.path.join(root_dir, "build", "oed")
    config.artifact_dir = os.path.join(root_dir, "lib")
    oed_objects = []
    for oed_source in oed_sources:
        oed_objects.append(config.fc(oed_source))
    liboed = config.ar(oed_objects, "liboed.a")

    # Build CInt (C interface) library
    config.source_dir = os.path.join(root_dir, "libcint")
    config.build_dir = os.path.join(root_dir, "build", "libcint")
    config.artifact_dir = os.path.join(root_dir, "lib")
    cint_objects = []
    for cint_source in cint_sources:
        cint_objects.append(config.cc(cint_source, include_dirs=[config.source_dir, os.path.join(root_dir, "include")]))
    libcint = config.ar(cint_objects, "libcint.a")

    # Build unit & perf tests
    config.source_dir = os.path.join(root_dir, "test")
    config.build_dir = os.path.join(root_dir, "build", "test")
    config.artifact_dir = os.path.join(root_dir, "test")
    include_dirs = [os.path.join(root_dir, "include"), os.path.join(root_dir, "test")]
    screening_object = config.cc("screening.c", include_dirs=include_dirs)
    config.ccld([config.cc("unit-test.c", include_dirs=include_dirs), screening_object], "unit-test",
        lib_dirs=[os.path.join(root_dir, "lib")],
        libs=["cint", "erd", "oed"],
        lib_deps=[liberd, liboed, libcint])
    config.ccld([config.cc("perf-test.c", include_dirs=include_dirs + [os.path.join(root_dir, "erd")]), screening_object], "perf-test",
        lib_dirs=[os.path.join(root_dir, "lib")],
        libs=["cint", "erd", "oed"],
        lib_deps=[liberd, liboed, libcint])


if __name__ == "__main__":
    sys.exit(main())
