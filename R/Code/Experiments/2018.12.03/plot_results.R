library(ggplot2)
library(reshape2)

load('/Users/tyler/RandomerForest/R/Results/2018.12.03/gmm_similarity.RData')

df <- rbind(data.frame(Alg=rep('SmerF', 3L), n.train=c(20L, 40L, 80L), MAE=apply(mae$SmerF, 2L, mean), SEM=apply(mae$SmerF, 2L, function(a) sd(a)/sqrt(length(a)))),
            data.frame(Alg=rep('RFD', 3L), n.train=c(20L, 40L, 80L), MAE=apply(mae$RFD, 2L, mean), SEM=apply(mae$RFD, 2L, function(a) sd(a)/sqrt(length(a)))))

p <- ggplot(data = df, aes(x=n.train, y=MAE, color=Alg)) +
  geom_line(size=2) +
  geom_errorbar(aes(ymin=MAE-SEM, ymax=MAE+SEM), size=2) +
  xlab("n.train") +
  ylab("MAE") +
  ggtitle("GMM") +
  scale_color_discrete() +
  guides(color = guide_legend(title = "Alg")) +
  theme_classic() +
  theme(text = element_text(size=20))

ggsave(filename='/Users/tyler/RandomerForest/R/Results/2018.12.03/gmm_mae.png', plot=p)

df <- rbind(data.frame(Alg=rep('SmerF', 3L), n.train=c(20L, 40L, 80L), train.time=apply(train.time$SmerF, 2L, mean), SEM=apply(train.time$SmerF, 2L, function(a) sd(a)/sqrt(length(a)))),
            data.frame(Alg=rep('RFD', 3L), n.train=c(20L, 40L, 80L), train.time=apply(train.time$RFD, 2L, mean), SEM=apply(train.time$RFD, 2L, function(a) sd(a)/sqrt(length(a)))))

p <- ggplot(data = df, aes(x=n.train, y=train.time, color=Alg)) +
  geom_line(size=2) +
  geom_errorbar(aes(ymin=train.time-SEM, ymax=train.time+SEM), size=2) +
  xlab("n.train") +
  ylab("train time (sec)") +
  ggtitle("GMM") +
  scale_color_discrete() +
  guides(color = guide_legend(title = "Alg")) +
  theme_classic() +
  theme(text = element_text(size=20))

ggsave(filename='/Users/tyler/RandomerForest/R/Results/2018.12.03/gmm_train_time.png', plot=p)

load('/Users/tyler/RandomerForest/R/Results/2018.12.03/radial_similarity.RData')

df <- rbind(data.frame(Alg=rep('SmerF', 3L), n.train=c(20L, 40L, 80L), MAE=apply(mae$SmerF, 2L, mean), SEM=apply(mae$SmerF, 2L, function(a) sd(a)/sqrt(length(a)))),
            data.frame(Alg=rep('RFD', 3L), n.train=c(20L, 40L, 80L), MAE=apply(mae$RFD, 2L, mean), SEM=apply(mae$RFD, 2L, function(a) sd(a)/sqrt(length(a)))))

p <- ggplot(data = df, aes(x=n.train, y=MAE, color=Alg)) +
  geom_line(size=2) +
  geom_errorbar(aes(ymin=MAE-SEM, ymax=MAE+SEM), size=2) +
  xlab("n.train") +
  ylab("MAE") +
  ggtitle("Radial Similarity") +
  scale_color_discrete() +
  guides(color = guide_legend(title = "Alg")) +
  theme_classic() +
  theme(text = element_text(size=20))

ggsave(filename='/Users/tyler/RandomerForest/R/Results/2018.12.03/radial_mae.png', plot=p)

df <- rbind(data.frame(Alg=rep('SmerF', 3L), n.train=c(20L, 40L, 80L), train.time=apply(train.time$SmerF, 2L, mean), SEM=apply(train.time$SmerF, 2L, function(a) sd(a)/sqrt(length(a)))),
            data.frame(Alg=rep('RFD', 3L), n.train=c(20L, 40L, 80L), train.time=apply(train.time$RFD, 2L, mean), SEM=apply(train.time$RFD, 2L, function(a) sd(a)/sqrt(length(a)))))

p <- ggplot(data = df, aes(x=n.train, y=train.time, color=Alg)) +
  geom_line(size=2) +
  geom_errorbar(aes(ymin=train.time-SEM, ymax=train.time+SEM), size=2) +
  xlab("n.train") +
  ylab("train time (sec)") +
  ggtitle("Radial Similarity") +
  scale_color_discrete() +
  guides(color = guide_legend(title = "Alg")) +
  theme_classic() +
  theme(text = element_text(size=20))

ggsave(filename='/Users/tyler/RandomerForest/R/Results/2018.12.03/radial_train_time.png', plot=p)


load('/Users/tyler/RandomerForest/R/Results/2018.12.03/angular_similarity.RData')

df <- rbind(data.frame(Alg=rep('SmerF', 3L), n.train=c(20L, 40L, 80L), MAE=apply(mae$SmerF, 2L, mean), SEM=apply(mae$SmerF, 2L, function(a) sd(a)/sqrt(length(a)))),
            data.frame(Alg=rep('RFD', 3L), n.train=c(20L, 40L, 80L), MAE=apply(mae$RFD, 2L, mean), SEM=apply(mae$RFD, 2L, function(a) sd(a)/sqrt(length(a)))))

p <- ggplot(data = df, aes(x=n.train, y=MAE, color=Alg)) +
  geom_line(size=2) +
  geom_errorbar(aes(ymin=MAE-SEM, ymax=MAE+SEM), size=2) +
  xlab("n.train") +
  ylab("MAE") +
  ggtitle("Angular Similarity") +
  scale_color_discrete() +
  guides(color = guide_legend(title = "Alg")) +
  theme_classic() +
  theme(text = element_text(size=20))

ggsave(filename='/Users/tyler/RandomerForest/R/Results/2018.12.03/angular_mae.png', plot=p)

df <- rbind(data.frame(Alg=rep('SmerF', 3L), n.train=c(20L, 40L, 80L), train.time=apply(train.time$SmerF, 2L, mean), SEM=apply(train.time$SmerF, 2L, function(a) sd(a)/sqrt(length(a)))),
            data.frame(Alg=rep('RFD', 3L), n.train=c(20L, 40L, 80L), train.time=apply(train.time$RFD, 2L, mean), SEM=apply(train.time$RFD, 2L, function(a) sd(a)/sqrt(length(a)))))

p <- ggplot(data = df, aes(x=n.train, y=train.time, color=Alg)) +
  geom_line(size=2) +
  geom_errorbar(aes(ymin=train.time-SEM, ymax=train.time+SEM), size=2) +
  xlab("n.train") +
  ylab("train time (sec)") +
  ggtitle("Angular Similarity") +
  scale_color_discrete() +
  guides(color = guide_legend(title = "Alg")) +
  theme_classic() +
  theme(text = element_text(size=20))

ggsave(filename='/Users/tyler/RandomerForest/R/Results/2018.12.03/angular_train_time.png', plot=p)
