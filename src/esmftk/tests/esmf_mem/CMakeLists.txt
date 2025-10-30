add_executable(esmf_mem mem_app.F90
                        mem_driver.F90
                        mem_comp1.F90
                        mem_comp2.F90
)
target_link_libraries(esmf_mem PUBLIC ESMF)
