
configure_file(mdcore_config.h.in mdcore_config.h @ONLY)

set (PUBLIC_HEADERS
  mdcore_f.h
  mdcore_d.h
  mdcore_config.h
  angle.h
  cutil_math.h
  errs.h
  mainpage.h
  potential_eval.h
  runner.h
  spme.h
  bond.h
  cycle.h
  exclusion.h
  queue.h
  runner_cuda.h
  task.h
  btree.h
  dihedral.h
  fptype.h
  part.h
  reader.h
  runner_cuda_main.h
  cell.h
  engine.h
  lock.h
  potential.h
  rigid.h
  space.h
  )

message("PUBLIC_HEADERS: ${PUBLIC_HEADERS}")

set(PUBLIC_HEADERS ${PUBLIC_HEADERS} PARENT_SCOPE)


