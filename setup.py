from setuptools import setup

with open('requirements.txt', 'r') as reqfile:
    requirements = [x.strip() for x in reqfile if x.strip()]

setup(
        name='chemworkflows',
        version='0.0.0-alpha1',
        packages=['chemworkflows', 'chemworkflows.apps'],
        url='',
        license='Apache 2.0',
        author='Aaron Virshup',
        author_email='aaron.virshup@autodesk.com',
        install_requires=requirements,
        description='',
        entry_points={
            'console_scripts': [
                'chemworkflow = chemworkflows.__main__:main'
            ]
        }
)
