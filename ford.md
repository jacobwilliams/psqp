project: psqp
src_dir: ./src
output_dir: ./doc
media_dir: ./media
project_github: https://github.com/jacobwilliams/psqp
summary: Modern Fortran edition of PSQP
author: Jacob Williams
github: https://github.com/jacobwilliams
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
         private
         protected
source: true
graph: true
search: true
preprocessor: gfortran -E
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html

{!README.md!}