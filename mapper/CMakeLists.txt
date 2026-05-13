cmake_minimum_required(VERSION 3.14)
project(mapper VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add bioparser
add_subdirectory(external/bioparser)

add_executable(${PROJECT_NAME}
    src/main.cpp
    src/mapper.cpp
    alignment/alignment.cpp
    minimizers/minimizers.cpp
)

# Link bioparser
target_link_libraries(${PROJECT_NAME}
    bioparser
)

# Include directories
target_include_directories(${PROJECT_NAME}
    PRIVATE external/bioparser/include
)

# Version definition for --version
target_compile_definitions(${PROJECT_NAME}
    PRIVATE PROJECT_VERSION="${PROJECT_VERSION}"
)