# Human_predict
Human navigational intent prediction

Although human navigational intent inference has been studied in the literature, none have adequately considered both the dynamics that describe human motion and internal human parameters that may affect human navigational behaviour.
A  general probabilistic framework to infer the probability distribution over future navigational states of a human. Our framework incorporates an extended Dubins car dynamics to model human movement, which captures differences in human navigational behaviour depending on their position, heading, and movement speed. We assume a noisily rational model of human behaviour that incorporates a) human navigational intent that may change over time, b) how optimal a person's actions are given the navigational intent, and c) how far ahead in time a person considers when choosing navigational actions.  These parameters are recursively and continuously updated in a Bayesian fashion. To make the Bayesian update and inference tractable, we exploit properties of the time-to-reach value function from optimal control and the extended Dubins car dynamics to construct a utility function on which the human policy is based, and employ particle representations of probability distributions where necessary.


How to refer:
::
@inproceedings{agandhuman,
  title={Human Navigational Intent Inference with Probabilistic and Optimal Approaches},
  author={Agand, Pedram and Taherahmadi, Mahdi and Lim, Angelica and Chen, Mo}
  booktitle={2022 IEEE International Conference on Robotics and Automation (ICRA)},
  pages={?},
  year={2022},
  organization={IEEE}
}
