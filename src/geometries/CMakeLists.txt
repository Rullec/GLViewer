set(
    geo_src
    Primitives.cpp
    Triangulator.cpp
)

add_library(geo_lib STATIC ${geo_src})