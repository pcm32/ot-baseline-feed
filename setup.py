from setuptools import setup, find_packages


def readme():
    with open('README.md') as f:
        return f.read()


print(find_packages())

setup(
        name='atlas-baseline-ot-feed',
        version='0.0.1',
        description='Expression Atlas baseline Open Targets exporter',
        long_description=readme(),
        packages=find_packages(),
        install_requires=['pandas'],
        author='Pablo Moreno',
        long_description_content_type='text/markdown',
        author_email='',
        scripts=['ot-export-atlas-baseline-exp.py'],
        license='MIT'
    )