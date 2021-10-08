|PyPI| |Docs|

.. |PyPI| image:: https://img.shields.io/pypi/v/tadmap.svg
   :target: https://pypi.org/project/tadmap
.. |Docs| image:: https://readthedocs.org/projects/tadmap/badge/?version=latest
   :target: https://tadmap.readthedocs.io/en/latest/?badge=latest



TAD Map - Consensus estimates of TAD layouts with associated tools for single-cell RNA-seq analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Topologically associating domains (TADs) are contiguous segments of the genome where the genomic elements are in frequent contact with each other. Genes that cooccupy a TAD are often functionally related and more likely to be coexpressed than random gene pairs, The TAD Map (preprint to come) provides a consensus estimate, aggregated from multiple experimental datasets, of this layout in human and mouse. It also provides tools to map any single-cell RNA-seq dataset to *TAD signatures*, where gene expression is mapped to TAD activation probabilities in each cell.

Read the documentation_.
We encourage you to report issues at our `Github page`_ ; you can also create pull reports there to contribute your enhancements.
If the TAD Map is useful in your research, please consider citing our preprint `bioRxiv (2021)`_.

.. _documentation: https://tadmap.readthedocs.io/en/latest/overview.html
.. _bioRxiv (2021): http://doi.org/10.1101/TBD
.. _Github page: https://github.com/rs239/tadmap
