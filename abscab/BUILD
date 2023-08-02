load("@rules_cc//cc:defs.bzl", "cc_library")

cc_library(
    name = "abscab",
    srcs = ["abscab.cc"],
    hdrs = ["abscab.hh"],
    visibility = ["//visibility:public"],
)

filegroup(
    name = "test_data",
    srcs = [
    	"resources/CircularWireLoop_A_phi_ref.dat",
    	"resources/idxRpStraightWireSegment.dat",
    	"resources/StraightWireSegment_B_phi_ref.dat",
    	"resources/testKnotsZpStraightWireSegment.dat",
    	"resources/testPointsStraightWireSegment.dat",
    	"resources/CircularWireLoop_B_rho_ref.dat",
    	"resources/idxZpCircularWireLoop.dat",
    	"resources/testKnotsRpCircularWireLoop.dat",
    	"resources/testPointsCircularWireLoop.dat",
    	"resources/testPointsZpCircularWireLoop.dat",
    	"resources/CircularWireLoop_B_z_ref.dat",
    	"resources/idxZpStraightWireSegment.dat",
    	"resources/testKnotsRpStraightWireSegment.dat",
    	"resources/testPointsRpCircularWireLoop.dat",
    	"resources/testPointsZpStraightWireSegment.dat",
    	"resources/idxRpCircularWireLoop.dat",
    	"resources/StraightWireSegment_A_z_ref.dat",
    	"resources/testKnotsZpCircularWireLoop.dat",
    	"resources/testPointsRpStraightWireSegment.dat",
    ],
)

cc_test(
    name = "TestAbscab",
    size = "small",
    srcs = [
    	"test_abscab.cc",
    ],
    visibility = ["//visibility:private"],
    data = [":test_data"],
    deps = [
    	":abscab",
    	"//abscab-cpp/test_util:test_util",
        "@com_google_googletest//:gtest_main",
        "@bazel_tools//tools/cpp/runfiles",
    ],
    copts = [
    	"-O3",
    	"-march=native",
    	"-fopenmp"
    ],
    linkopts = [
        "-lgomp",
    ],
)

cc_binary(
	name = "DemoAbscab",
	srcs = [
		"demo_abscab.cc"
	],
	visibility = ["//visibility:private"],
	deps = [
    	":abscab",
    	"//abscab-cpp/test_util:test_util",
    ],
    copts = [
    	"-O3",
    	"-march=native",
    	"-fopenmp"
    ],
    linkopts = [
        "-lgomp",
    ],
)
