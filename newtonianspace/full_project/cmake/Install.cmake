include(GNUInstallDirs)
# NewtonianSpace library
install(TARGETS newtonianspace_obj newtonianspace_shared newtonianspace_static
  EXPORT NewtonianSpaceLibrary
  ARCHIVE COMPONENT development
  LIBRARY COMPONENT runtime
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/newtonianspace
    COMPONENT runtime
)

if (UNIX)
  install(CODE "execute_process(COMMAND ldconfig)"
          COMPONENT runtime
  )
endif()

install(EXPORT NewtonianSpaceLibrary
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/newtonianspace/cmake
  #NAMESPACE NewtonianSpace::
  COMPONENT runtime
)

install(FILES "NewtonianSpaceConfig.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/newtonianspace/cmake
)

# CalcConsole runtime
#[[install(TARGETS calc_console
  RUNTIME COMPONENT runtime
)]]

# CPack configuration
set(CPACK_PACKAGE_VENDOR "Rabometrics")
set(CPACK_PACKAGE_CONTACT "rabometrics@gmail.com")
set(CPACK_PACKAGE_DESCRIPTION "Newtonian Space")
include(CPack)
