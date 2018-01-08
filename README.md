
  What is it?
  -----------
  The zip folder should contain three python scripts, 'basicSGZ.py' is a basic method to predict germline/somatic
  for short variants; 'fmiSGZ.py' is the method developed in Foundation Medicine Inc to predict germline/somatic origins
  and zygosity for short variants; 'run_test.py' is a script to test the previous two scripts on a test sample set.

  Installation and how to run
  ------------------------

  FMI SGZ method does not require installation. The core method is implemented in python
  script fmiSGZ.py. Simply provide required files while calling the method in python terminals.

  After download and extract all files, run script 'run_test.py' to test if the scripts are working
  as expected. 'run_test.py' runs the basic SGZ method and FMI SGZ method on four test samples. If the process
  succeeded, a message '-------Test succeeded.-------' should be printed on the standard output.

  All the scripts are developed under Python 2.7.6.

  Contacts
  --------

  Please contact James Sun (jsun@foundationmedicine.com) if you have any questions.
