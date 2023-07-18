load("@rules_cc//cc:defs.bzl", "cc_library")

cc_library(
    name = "cel",
    srcs = ["cel.cc"],
    hdrs = ["cel.hh"],
    visibility = ["//:__subpackages__"],
    copts = [
    	"-O3",
    	"-march=native",
    	"-fopenmp",
    ],
)

cc_test(
    name = "TestCel",
    size = "small",
    srcs = [
    	"test_cel.cc",
    ],
    visibility = ["//visibility:private"],
    deps = [
    	":cel",
    	"//abscab-cpp/util:util",
        "@com_google_googletest//:gtest_main",
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

