# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "C:/Users/maygb/source/repos/mboffelli/klue/ext/bifrost"
  "C:/Users/maygb/source/repos/mboffelli/klue/ext/bifrost/src/bifrost-build"
  "C:/Users/maygb/source/repos/mboffelli/klue/ext/bifrost"
  "C:/Users/maygb/source/repos/mboffelli/klue/ext/bifrost/tmp"
  "C:/Users/maygb/source/repos/mboffelli/klue/ext/bifrost/src/bifrost-stamp"
  "C:/Users/maygb/source/repos/mboffelli/klue/ext/bifrost/src"
  "C:/Users/maygb/source/repos/mboffelli/klue/ext/bifrost/src/bifrost-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "C:/Users/maygb/source/repos/mboffelli/klue/ext/bifrost/src/bifrost-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "C:/Users/maygb/source/repos/mboffelli/klue/ext/bifrost/src/bifrost-stamp${cfgdir}") # cfgdir has leading slash
endif()
