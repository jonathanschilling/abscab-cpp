load("@rules_cc//cc:defs.bzl", "cc_library")

cc_library(
    name = "compsum",
    srcs = ["compsum.cc"],
    hdrs = ["compsum.h"],
    visibility = ["//:__subpackages__"],
    copts = [
    	"-O3",
    	"-march=native",
    	"-fopenmp",
    ],
)
