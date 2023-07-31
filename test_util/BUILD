load("@rules_cc//cc:defs.bzl", "cc_library")

cc_library(
    name = "test_util",
    srcs = ["test_util.cc"],
    hdrs = ["test_util.hh"],
    visibility = ["//:__subpackages__"],
    copts = [
    	"-O3",
    	"-march=native",
    	"-fopenmp"
    ],
)
