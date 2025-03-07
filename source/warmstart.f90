subroutine warmstart ()

  use parameters
  use globals

  implicit none

  character(len = 200)   :: fp

  write(fp,'((a),i3.3)') '/storage5/luis.arcos/out200SG/200',rank
