all: SampleArcLenCF.xml CurvatureStreamlinesCF.xml JoinStreamlinesCF.xml OrientedGlyphCF.xml SampleSurfaceCF.xml MetropolisHastingsCF.xml  SampleStreamlinesInTime.xml

SampleStreamlinesInTime.xml: SampleStreamlinesInTime.py
	python python_filter_generator.py SampleStreamlinesInTime.py SampleStreamlinesInTime.xml

MetropolisHastingsCF.xml: MetropolisHastingsCF.py
	python python_filter_generator.py MetropolisHastingsCF.py MetropolisHastingsCF.xml

OrientedGlyphCF.xml: OrientedGlyphCF.py
	python python_filter_generator.py OrientedGlyphCF.py OrientedGlyphCF.xml

JoinStreamlinesCF.xml: JoinStreamlinesCF.py
	python python_filter_generator.py JoinStreamlinesCF.py JoinStreamlinesCF.xml

SampleArcLenCF.xml: SampleArcLenCF.py
	python python_filter_generator.py SampleArcLenCF.py SampleArcLenCF.xml

CurvatureStreamlinesCF.xml: CurvatureStreamlinesCF.py
	python python_filter_generator.py CurvatureStreamlinesCF.py CurvatureStreamlinesCF.xml

SampleSurfaceCF.xml: SampleSurfaceCF.py
	python python_filter_generator.py SampleSurfaceCF.py SampleSurfaceCF.xml
