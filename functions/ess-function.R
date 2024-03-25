# Efective Sample Size Function
# inputs: data frame, treatment column in "", and weights column in ""

ess <- function(data, trt, wts) {
  data.trt <- data %>% filter(.data[[trt]] == 1)
  n.trt <-(sum(data.trt[[wts]])^2)/ (sum(data.trt[[wts]]^2))
  data.ctrl <- data %>% filter(.data[[trt]] == 0)
  n.ctrl <- (sum(data.ctrl[[wts]])^2)/ (sum(data.ctrl[[wts]]^2))
  n.total <- n.trt + n.ctrl
  
  result <- data.frame(
    "N treated" = n.trt,
    "N control" = n.ctrl,
    "N total" = n.total
  )
  
  return(result)
}




