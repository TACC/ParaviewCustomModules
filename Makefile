
## ========================================================================== ##
## Copyright (c) 2019 The University of Texas at Austin.                      ##
## All rights reserved.                                                       ##
##                                                                            ##
## Licensed under the Apache License, Version 2.0 (the "License");            ##
## you may not use this file except in compliance with the License.           ##
## A copy of the License is included with this software in the file LICENSE.  ##
## If your copy does not contain the License, you may obtain a copy of the    ##
## License at:                                                                ##
##                                                                            ##
##     https://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
## Unless required by applicable law or agreed to in writing, software        ##
## distributed under the License is distributed on an "AS IS" BASIS, WITHOUT  ##
## WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.           ##
## See the License for the specific language governing permissions and        ##
## limitations under the License.                                             ##
##                                                                            ##
## ========================================================================== ##

all: SampleArcLenCF.xml CurvatureStreamlinesCF.xml JoinStreamlinesCF.xml OrientedGlyphCF.xml SampleSurfaceCF.xml MetropolisHastingsCF.xml  SampleStreamlinesInTime.xml TrimStreamlinesCF.xml MH.xml MH2.xml MH3.xml DecimateStreamlines.xml MH4.xml PDF.xml

DecimateStreamlines.xml: DecimateStreamlines.py
	python2 python_filter_generator.py DecimateStreamlines.py DecimateStreamlines.xml

PDF.xml: PDF.py
	python2 python_filter_generator.py PDF.py PDF.xml

TrimStreamlinesCF.xml: TrimStreamlinesCF.py
	python2 python_filter_generator.py TrimStreamlinesCF.py TrimStreamlinesCF.xml

SampleStreamlinesInTime.xml: SampleStreamlinesInTime.py
	python2 python_filter_generator.py SampleStreamlinesInTime.py SampleStreamlinesInTime.xml

MetropolisHastingsCF.xml: MetropolisHastingsCF.py
	python2 python_filter_generator.py MetropolisHastingsCF.py MetropolisHastingsCF.xml

MH.xml: MH.py
	python2 python_filter_generator.py MH.py MH.xml

MH2.xml: MH2.py
	python2 python_filter_generator.py MH2.py MH2.xml

MH3.xml: MH3.py
	python2 python_filter_generator.py MH3.py MH3.xml

MH4.xml: MH4.py
	python2 python_filter_generator.py MH4.py MH4.xml

OrientedGlyphCF.xml: OrientedGlyphCF.py
	python2 python_filter_generator.py OrientedGlyphCF.py OrientedGlyphCF.xml

JoinStreamlinesCF.xml: JoinStreamlinesCF.py
	python2 python_filter_generator.py JoinStreamlinesCF.py JoinStreamlinesCF.xml

SampleArcLenCF.xml: SampleArcLenCF.py
	python2 python_filter_generator.py SampleArcLenCF.py SampleArcLenCF.xml

CurvatureStreamlinesCF.xml: CurvatureStreamlinesCF.py
	python2 python_filter_generator.py CurvatureStreamlinesCF.py CurvatureStreamlinesCF.xml

SampleSurfaceCF.xml: SampleSurfaceCF.py
	python2 python_filter_generator.py SampleSurfaceCF.py SampleSurfaceCF.xml
