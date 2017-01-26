from distutils.core import setup

setup(
        name='chemworkflows',
        version='0.0.0-alpha1',
        packages=['chemworkflows', 'chemworkflows.apps'],
        url='',
        license='Apache 2.0',
        author='Aaron Virshup',
        author_email='aaron.virshup@autodesk.com',
        description='',
        entry_points={
            'console_scripts': [
                'chemworkflow = chemworkflows.__main__:main'
            ]
        }
)
