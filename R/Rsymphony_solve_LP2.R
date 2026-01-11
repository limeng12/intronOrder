# 最简化的智能选择方案
Rsymphony_solve_LP2 <- function(obj, mat, dir, rhs, types = NULL, max = FALSE) {
  
  # 尝试Rsymphony
  if (requireNamespace("Rsymphony", quietly = TRUE) ) {
    #print("using Rsymphony")
    return(Rsymphony::Rsymphony_solve_LP(obj=obj, mat=mat, dir=dir, rhs=rhs, 
                                         types=types, max=max) );
    
  }
  
  # 尝试lpSolve
  if (requireNamespace("lpSolve", quietly = TRUE) ) {
    #print("using lpSolve")
    
    return(.simple_lpSolve_wrapper(obj, mat, dir, rhs, types, max) );
    
  }
  
  # 安装lpSolve
  # install.packages("lpSolve", repos = "https://cloud.r-project.org", quiet = TRUE)
  # library(lpSolve, quietly = TRUE)
  # return(.simple_lpSolve_wrapper(obj, mat, dir, rhs, types, max))
}

.simple_lpSolve_wrapper <- function(obj, mat, dir, rhs, types, max) {
  result <- lpSolve::lp(
    direction = ifelse(max, "max", "min"),
    objective.in = obj,
    const.mat = mat,
    const.dir = dir,
    const.rhs = rhs,
    all.bin = if (!is.null(types) && all(types == "B")) TRUE else FALSE
  )
  
  return(list(
    objval = result$objval,
    solution = result$solution,
    status = result$status
  ))
}


