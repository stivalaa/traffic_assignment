from distutils.core import setup, Extension
 
module1 = Extension('cdijkstra',
                    include_dirs = ['libpqueue/src'],
                    sources = ['dijkstra_module.c',
                               'libpqueue/src/pqueue.c'])

module2 = Extension('pape',
                    sources = ['pape_module.c'])

setup (name = 'Dijkstra',
        version = '1.0',
        description = "Shortest paths using Dijkstra's algorithm",
        ext_modules = [module1,module2])
