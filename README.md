# dNNDR

> Deep neural network assisted drug recommendation (dNNDR) system

Introduction
------------

Identification of novel drug target interactions is a labor intensive, low throughput process. In-silico alternatives have proved to be of immense importance in assisting the drug discovery process. In this paper, we present dNNDR, an att-biLSTM and gCNN based multi-label classification system that segregates drug-target pairs into three categories. SMILEs and proteins sequences, along with molecular descriptors are embedded into a machine interpretable form to extract critical features form the ordered data.  We used KIBA and Davis datasets for training and external validation respectively. Tools such as CHEMBL web resource, iFeature and in-house developed dNNDR-featx were employed for data retrieval and processing. Experimental results suggest an improvement in performance in terms of auROC, auPR, F1-score and accuracy over baseline methods in internal as well as external validation sets over baseline methods. Making use of the ordered information present in datasets like SMILEs and protein sequence was central to the idea of dNNDR. Looking at the exceptional performance on an external dataset suggest that the proposed methodology is able to make sense of the data and elucidate critical features that govern the relationship between a drug and its target. Additionally, we also introduce some novel drug-target pairs form gold standard enzyme dataset demonstrating a practical use case.


Naming convention: **X DXXX TXX** | **12 3456 78**

1. C: *Classification*, R: *Regression*
2. D: *Drug*
3. Sequence *(Included-1, Excluded-0)*
4. XAE *(Included-1, Excluded-0)*
5. Descriptors *(Included-1, Excluded-0)*
6. T: *Target*
7. Sequence *(Included-1, Excluded-0)*
8. Descriptors *(Included-1, Excluded-0)*

Example: **CD100T010**

Data/Notebook for a classification model with Drug features containing only **XAE** and target features having only **Descriptors** as features.

Usage
------------

* **Modular approach** — We propose to harness the potential of state-of-the-art ML (Machine Learning) algorithms like Bi-directional Long Short-Term Memory (bi-LSTM) and Graph Convolutional neural networks (gCNNs)

Questions? Need Help? Found a bug?
--------------------
Found a bug with upstream Slate? Go ahead and [submit an issue](https://github.com/iamysk/dNNDR/issues). And, of course, feel free to submit pull requests with bug fixes or changes to the `dev` branch.

Contributors
--------------------
- [Shashank Yadav](https://github.com/xinformatics)
- [Yogesh Kalakoti](https://github.com/iamysk)
