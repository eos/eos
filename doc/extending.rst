#############
Extending EOS
#############

The three most common cases of extending are:

* adding a new parameter;
* adding a new measurement or theory prediction as a constraint;
* adding a new theory expression for a (pseudo) observable.

To make your modification to EOS permanent, we ask you to submit
it as a Github pull request.

**********************
Adding a new parameter
**********************

Parameters are stored in YAML files below eos/parameters. Each parameter entry consists of the
parameter's name, its default value and range, as well as optionally further metadata.

To be able to create a summary of all EOS parameters on this website, the format of the
parameter files is hierarchical.

Each parameter file corresponds to a section, which is
subdivided into groups. Each section is a key/value map with the following mandatory keys:

* **title** (*str*) -- The section's title. Basic LaTeX math mode elements are permitted, provided that
  they are used within inline mathmode (i.e., enclosed by single dollar signs).
* **description** (*str*) -- The section's (succinct) description. Basic LaTeX math mode elements are permitted as for
  a section title.
* **groups** (*list*) -- The list of parameter groups in this section.

Each parameter group's entry is a key/value map with the following mandatory keys:

* **title** (*str*) -- The group's title.  Basic LaTeX math mode elements are permitted as for
  a section title.
* **description** (*str*) -- The group's (succinct) descripton. Basic LaTeX math mode elements are permitted as for
  a section title.
* **parameters** (*map*) -- The group's parameters. The keys into this map are strings representing a valid
  :class:`QualifiedName <eos.QualifiedName>`. The values are parameter entries.
  
Each parameter' entry is a key/value map with the following mandatory keys:

* **central** (*float*) -- The default central value of the parameter.
* **min** (*float*) -- The minimal value the parameter can take.
* **max** (*float*) -- The maximal value the parameter can take.

Further optional but recommended keys are:

* **latex** (*str*) -- The LaTeX representation of the parameter in math mode.
  Enclosing math mode delimiters (such as ``$...$``, ``$$...$$``, or ``\(...\)``) are not permitted.
* **unit** (*str*) -- The unit in which the parameter value is expressed.
* **comment** (*str*) -- A comment on the source of the parameter's default value.

New parameters should be added to existing groups best as possible.

Example
#######

.. code-block::

   title: Example parameter section for the documentation
   description: This section collects parameters of groups A and B
   groups:
     - title: Example group A
       description: Group A's parameter are relevant to a a specific process $X\to Y$, and grouped together.
       parameters:
         X->Y::a0@B:2002A
           central: +1.0
           min:      0.0
           max:     +2.0
           latex:   'a_0^{X\to Y}'
           # arbitrary units
           comment: 'Taken from ref. [B:2002A] without changes.'
        X->Y::a2@B:2002A
           central: -0.1
           min:     -2.0
           max:     +2.0
           latex:   'a_1^{X\to Y}'
           # arbitrary units
           comment: 'Taken from ref. [B:2002A] without changes.'
     - title: Example group B
       description: Group B's parameter are also relevant to a a specific process $X\to Y$, but are independent of group A
       parameters:
         X->Y::rho
           central: +1.234
           min:     -2.000
           max:     +2.000
           latex:    '\rho(X\to Y)'
           comment:  'relevant only in certain corner cases.'


