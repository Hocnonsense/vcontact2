from setuptools import setup

setup(name='vcontact2',
      version='0.9.2b',
      description='Viral Contig Automatic Clutering and Taxonomy',
      url='https://bitbucket.org/MAVERIClab/vcontact2',
      author='Benjamin Bolduc',
      long_description_markdown_filename='README.md',
      author_email='bolduc.10@osu.edu',
      license='GPLv3',
      packages=['vcontact', 'vcontact.exports', 'vcontact.utilities'],
      package_data={'vcontact': ['data/ViralRefSeq-prokaryotes-v88.faa.gz',
                                 'data/ViralRefSeq-prokaryotes-v88.protein2contig.csv'
                                 'data/ViralRefSeq-prokaryotes-v88.Merged-reference.csv'
                                 'data/ViralRefSeq-prokaryotes-v85.faa.gz',
                                 'data/ViralRefSeq-prokaryotes-v85.protein2contig.csv',
                                 'data/ViralRefSeq-prokaryotes-v85.ICTV-reference.csv',
                                 'data/ViralRefSeq-prokaryotes-v85.Merged-reference.csv',
                                 'data/ViralRefSeq-archaea-v85.faa.gz',
                                 'data/ViralRefSeq-archaea-v85.protein2contig.csv',
                                 'data/ViralRefSeq-archaea-v85.Merged-reference.csv'
                                 ]},
      scripts=['bin/vcontact'],
      setup_requires=['setuptools-markdown'],
      install_requires=[
        'networkx>=1.11',
        'numpy>=1.12.1',
        'scipy>=0.19.0',
        'pandas>=0.21.0',
        'scikit-learn>=0.18.1',
        'biopython>=1.68',
        'hdf5>=1.8.17',
        'tables>=3.3.0'
      ]
      )
