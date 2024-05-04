The Data Model
==============

Reaction
--------
A reaction object captures a single chemical or biological reaction.

.. py:class:: rbc2.Reaction

   A dataclass representing a reaction.  Reactions can be combined into Pathways

   .. py:attribute:: product

      A string representing the product of the reaction.

   .. py:attribute:: substrates

      A list of strings representing the substrates of the reaction.

   .. py:attribute:: unique_id

      A unique identifier for the reaction. If not provided, a UUID will be generated.

   .. py:attribute:: name

      The name of the reaction.

   .. py:attribute:: rxn_type

      The type of the reaction.

   .. py:attribute:: rxn_domain

      The domain of the reaction.  This will typically be 'biocatalysis', 'biosynthesis' or 'chemistry'

   .. py:attribute:: score

      A score for the reaction. Defaults to 0.

   .. py:attribute:: template_metadata

      A dictionary with metadata for the reaction template. Defaults to an empty dictionary.

   .. py:attribute:: precedents

      A list of Precedent objects representing the precedents of the reaction. Defaults to an empty list.

   .. py:attribute:: feasability_filter_scores

      A dictionary of feasibility filter scores. Defaults to an empty dictionary.


   .. py:method:: get_complexity_change()

      Returns the change in complexity of the reaction. If not already calculated, it is calculated and stored.

   .. py:method:: get_similarity_score()

      Returns the best similarity score of the reaction from its precedents. If there are no precedents, it returns 0.

   .. py:method:: reaction_smiles()

      Returns a SMILES string representing the reaction.  For example 'CCO>>CC=O'

   .. py:method:: to_dict()

      Returns a dictionary representation of the Reaction object.