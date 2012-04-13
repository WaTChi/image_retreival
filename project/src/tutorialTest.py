# context and query objects we override
from context import _Context, _Query

# os utils
import os

class TestContext(_Context):
    """The context object contains all the parameters the query system
       takes, as well as pre-defined information about query sets.
       
       It is important to note that the Context object was originally used to
       make running experiments less painful and has many internal
       parameters that are used to drive experiments. We don't use many
       of these here - see context.py for more details.

       In this simple example we only care to override enough properties
       to make the pipeline use our test database."""

    def __init__(self):
        # Construct the parent class first
        super(TestContext, self).__init__()

        # Used internally by context to differentiate between query sets
        # Set to a dummy value for this example.
        self.QUERY = 'SingleImageTest'

    @property
    def dbdump(self):
        """Where the preprocessed database files are"""
        return 'testdb'

    @property
    def dbdir(self):
        """Where the cell configurations are"""
        return 'testdb/cells-236.6'

    @property
    def matchdir(self):
        """Temporary directory used for caching data"""
        p = 'test_matches'
        if not os.path.exists(p):
            os.mkdir(p)
        return p

# utilities for extracting features from the query image
import client

# framework for running queries
import system

# This is the query image
image = 'tutorial/example_query_images/query1.png'

# We need to extract features from the incoming query image
client.preprocess_image(image, image[:-4] + '.pgm', width=200, height=200)
client.extract_features(image[:-4] + '.pgm')

# Construct a context
C = TestContext()
C.check()

# Construct query object
Q = _Query()
Q.jpgpath = image
Q.siftpath = os.path.splitext(image)[0] + "sift.txt"
Q.setSensorCoord(37.875507, -122.264883)
Q.check()

# Tell the system to match this query Q in the context C
stats, matchedimg, matches, ranked = system.match(C, Q)

# The identifier for the matched database image
print "Matched db image ", matchedimg

# The system also spits out a db-query feature match visualization
print "Visualization in ", C.resultsdir
