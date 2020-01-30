library("lattice")

logistic <- function(x) {
    return(1.0 / (1.0 + exp(-x)))
}

# we have h(t+1) = h(t) * (1.0 - g) + u
# hence, h_eq = u/g

l_influx <- -1.5
l_feedback <- -1.0
l_stress_influx <- -4.0

print(logistic(l_influx) / logistic(l_feedback) + logistic(l_stress_influx))
