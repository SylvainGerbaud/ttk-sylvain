ttk_add_base_template_library(localizedTopologicalSimplificationMPI
  HEADERS
    LocalizedTopologicalSimplificationMPI.h
  DEPENDS
    triangulation
    localizedTopologicalSimplification
    Boost::boost
)