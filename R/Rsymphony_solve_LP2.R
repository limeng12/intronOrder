# 最简化的智能选择方案（优先 highs）
Rsymphony_solve_LP2 <- function(obj, mat, dir, rhs, types = NULL, max = FALSE) {
  
  # 1. 优先尝试 highs
  if (requireNamespace("highs", quietly = TRUE)) {
    # 转换类型格式：highs 接受 "C" / "I" / "B" 等
    highs_types <- NULL
    if (!is.null(types)) {
      highs_types <- ifelse(types == "B", "I",
                            ifelse(types == "I", "I",
                                   ifelse(types == "C", "C", "C")))
      # 二进制变量在 highs 中也是整数加边界 [0,1]
    }
    
    result <- highs::highs_solve(
      L = obj,
      lower = rep(0, length(obj)),   # 默认下界（可改进）
      upper = rep(Inf, length(obj)), # 默认上界
      A = mat,
      lhs = ifelse(dir == ">=", rhs, -Inf),
      rhs = ifelse(dir == "<=", rhs, Inf),
      types = highs_types,
      maximum = max
    )
    
    # 转换为统一输出格式
    return(list(
      objval = result$objective_value,
      solution = result$primal_solution,
      status = result$status
    ))
  }
  
  # 2. 其次尝试 Rsymphony
  if (requireNamespace("Rsymphony", quietly = TRUE)) {
    return(Rsymphony::Rsymphony_solve_LP(
      obj = obj, mat = mat, dir = dir, rhs = rhs,
      types = types, max = max
    ))
  }
  
  # 3. 最后尝试 lpSolve
  if (requireNamespace("lpSolve", quietly = TRUE)) {
    return(.simple_lpSolve_wrapper(obj, mat, dir, rhs, types, max))
  }
  
  stop("No solver available (highs, Rsymphony, lpSolve)")
}

# lpSolve 包装器（保持不变）
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



