

zermelo = function(game, win, lambda = 0.1, K = 300, tol = 0.001, beta = 0.0001, x0 = rep(1 / nrow(game), nrow(game))) {
  if (!is.matrix(game) & ! isSymmetric(game)){
    stop("'game' must be a symmetric matrix.")
  }
  
  n = nrow(game)
  
  if (!is.numeric(win) & length(win) != n) {
    stop("'win' must be a numeric vector of the same length as the row dimension of 'game'")
  }
  
  x = x0  
  y = numeric(n)
  k = 0; err = 1
  
  while(k < K & err > tol) {
    err = 0
    for (i in 1:n) {
      d = 0
      for (j in 1:n) {
        d = d + (game[i,j] + lambda)/ (x[i] + x[j] + beta)
      }
      y[i] = (win[i] + lambda*(n-1)/2) / d
    }
    y = y / sum(y)  
    err = sum(abs(y - x))
    x = y
    k = k + 1
  }  
  print(k)
  print(err)
  return(x)
}

zermelo.system = function(x) {
  teams = sort(unique(c(x$Winner, x$Loser)))
  n = length(teams)
  #initial game matrix
  game = matrix(0, nrow = n, ncol = n)
  #initial win vector
  win = numeric(n)
  w = ifelse(x$PD > 0 , 1, 0.5)
  
  for(k in 1:nrow(x)){
    i = match(x$Winner[k], teams)
    j = match(x$Loser[k], teams)
    
    game[i,j] = game[i,j] + 1
    game[j,i] = game[i,j]
    
    win[i] = win[i] + w[k]
    win[j] = win[j] + (1 - w[k])
  }
  
  return(list(game = game, win = win))
}


accuracy = function(r, game, win, type = 1, K1 = 0.00001, K2 = 0.001) {
  n = length(win)
  p = matrix(0, nrow = n, ncol = n)
  
  if (type == 1) {
    for (i in 1:n) {
      for (j in 1:n) {
        p[i,j] = r[i] / (r[i] + r[j] + K1)
      }
    }
  } else {
    for (i in 1:n) {
      for (j in 1:n) {
        p[i,j] = 1  / (1 + exp(-K2 * (r[i] - r[j])))
      }
    }
  }
  diag(p) = 0
  
  exp_win = numeric(n)
  for (i in 1:n) {
    exp_win[i] = sum(game[i,] * p[i,])  
  }
  sqrt(mean((win - exp_win)^2))
}









