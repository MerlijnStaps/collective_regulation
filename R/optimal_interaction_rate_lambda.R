source("R/plotting.R")
source("R/functions.R")



par(pty="m")
d1 = read.table("results/optimal_interaction_rates_lambda.csv",header=T,stringsAsFactors=F,sep=",")
d1 = add_random_switching(d1)

b0 = 2
n0 = 0
a0 = 1


strategies = c("individual_regulation","always_active","always_inactive","collective_regulation_I")

d = subset(d1,d1$b == b0 & d1$alpha == a0 & d1$n == n0)
d = add_random_switching(d, add_random_inactivation = FALSE)
d = calculate_winners(d, "collective_regulation_I")




pdf("plots/optimal_interaction_rate_lambda.pdf",width=4.5,height=4.5)
par(mgp=c(2,0.75,0))
interaction_rate_plot_lambda(d, analytical_predictions=T, background_colors = NA,
                            main = "Optimal interaction rate")
interaction_rate_plot_lambda(d, analytical_predictions=T, background_colors = strategy_colors[d$winner], 
                            main = "Optimal strategy")
interaction_rate_plot_lambda(d, analytical_predictions=T, background_colors = strategy_colors[d$random_switching_winner], 
                            main = "Optimal strategy\n(with spontaneous switching)")
interaction_rate_plot_lambda(d, analytical_predictions=F, background_colors = sapply(d$fitness_increase, function(x){fitness_increase_color(x, max_value = 0.08)}),
                            main = "Fitness increase due to interactions")
interaction_rate_plot_lambda(d, analytical_predictions=F, background_colors = sapply(d$random_switching_fitness_increase, function(x){fitness_increase_color(x, max_value = 0.08)}),
                            main = "Fitness increase due to interactions\n(with spontaneous switching)")
dev.off()

