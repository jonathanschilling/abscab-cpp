load("@rules_cc//cc:defs.bzl", "cc_library")

cc_library(
    name = "util",
    srcs = ["util.cc"],
    hdrs = ["util.h"],
    visibility = ["//:__subpackages__"],
    copts = [
    	"-O3",
    	"-march=native",
    	"-fopenmp"
    ],
)
