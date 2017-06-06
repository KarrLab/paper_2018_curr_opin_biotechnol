""" Analysis of the content of BioModels

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-31
:Copyright: 2017, Karr Lab
:License: MIT
"""

import bioservices
import ete3
import glob
import kinetic_datanator.data_source.bio_portal
import libsbml
import os
import sqlalchemy
import sqlalchemy.ext.declarative
import sqlalchemy.orm
import subprocess
import wc_utils.util.list
import wc_utils.workbook.core
import wc_utils.workbook.io

RELEASE_NAME = 'BioModels_Database-r30_pub-sbml_files'
SBML_FILES_URL = 'ftp://ftp.ebi.ac.uk/pub/databases/biomodels/releases/2016-05-10/{}.tar.bz2'.format(RELEASE_NAME)
DATA_DIRNAME = os.path.join(os.path.dirname(__file__), 'data')
SBML_FILES_ARCHIVE_FILENAME = os.path.join(DATA_DIRNAME, 'sbml_files.tar.bz2')
SBML_FILES_DIRNAME = os.path.join(DATA_DIRNAME, 'sbml_files')
SBML_FILES_DATABASE_FILENAME = os.path.join(DATA_DIRNAME, 'sbml_files.sqlite')
ANNOTATIONS_EXCEL_FILENAME = os.path.join(DATA_DIRNAME, 'models.xlsx')
SUMMARY_EXCEL_FILENAME = os.path.join(DATA_DIRNAME, 'summary.xlsx')


def create_data_directory():
    """ Create directory for files """
    if not os.path.isdir(DATA_DIRNAME):
        os.makedirs(DATA_DIRNAME)


def download_biomodels():
    """ Download BioModels release and extract content """

    # download release
    if not os.path.isfile(SBML_FILES_ARCHIVE_FILENAME):
        print('Downloading BioModels ...')
        subprocess.call(['wget', SBML_FILES_URL, '-O', SBML_FILES_ARCHIVE_FILENAME])
        print('  done.')

    # extract release archive
    if not os.path.isdir(SBML_FILES_DIRNAME):
        print('Unpacking BioModels ...')
        subprocess.call(['tar', '-xvjf', SBML_FILES_ARCHIVE_FILENAME, '-C', DATA_DIRNAME])
        os.rename(os.path.join(DATA_DIRNAME, RELEASE_NAME), SBML_FILES_DIRNAME)
        print('  done.')


def get_database_engine():
    """
    Returns:
        :obj:`sqlalchemy.engine.Engine`: database engine
    """
    return sqlalchemy.create_engine('sqlite:///' + SBML_FILES_DATABASE_FILENAME)


def get_database_session():
    """
    Returns:
        :obj:`sqlalchemy.orm.session.Session`: sqlalchemy session
    """
    engine = get_database_engine()
    return sqlalchemy.orm.sessionmaker(bind=engine)()


def setup_database(clear=False):
    if not os.path.isfile(SBML_FILES_DATABASE_FILENAME) or clear:
        clear_database()


def clear_database():
    engine = get_database_engine()
    Base.metadata.drop_all(engine)
    Base.metadata.create_all(engine)


def load_database():
    sbml_reader = libsbml.SBMLReader()
    session = get_database_session()

    if session.query(Model).count() > 0:
        return

    print('Loading models ...')

    curated_models = sorted(glob.glob(os.path.join(SBML_FILES_DIRNAME, 'curated', '*.xml')))
    for i_model, filename in enumerate(curated_models):
        if i_model % 100 == 0:
            print('  Loading curated model {} of {}'.format(i_model + 1, len(curated_models)))
        model = load_model_into_database(filename, True, sbml_reader, session)

    non_curated_models = sorted(glob.glob(os.path.join(SBML_FILES_DIRNAME, 'non_curated', '*.xml')))
    for i_model, filename in enumerate(non_curated_models):
        if i_model % 100 == 0:
            print('  Loading non-curated model {} of {}'.format(i_model + 1, len(non_curated_models)))
        model = load_model_into_database(filename, False, sbml_reader, session)

    print('  done.')

    session.commit()


def load_model_into_database(filename, curated, sbml_reader, session):
    """
    Args:
        filename (:obj:`str`): path to a SBML file
        curated (:obj:`bool`): :obj:`True`, the model has been curated
        sbml_reader (:obj:`libsbml.SBMLReader`): SBML file reader
        session (:obj:`sqlalchemy.orm.session.Session`): sqlalchemy session

    Returns:
        :obj:`Model`: model
    """
    # todo: detect mathematical type (ODE, SSA, logical, FBA, spatial, rule-based)

    doc = sbml_reader.readSBMLFromFile(filename)
    sbml_model = doc.getModel()
    if not sbml_model:
        return None

    id, _, _ = os.path.basename(filename).partition('.xml')
    label = sbml_model.getId()
    name = sbml_model.getName()
    annotations = parse_model_annotations(sbml_model, session)
    type = parse_model_type(doc, annotations)

    num_reaction_parameters = 0
    reactions_sbml = sbml_model.getListOfReactions()
    for i_reaction in range(sbml_model.getNumReactions()):
        reaction_sbml = reactions_sbml.get(i_reaction)
        kinetic_law_sbml = reaction_sbml.getKineticLaw()
        if kinetic_law_sbml:
            num_reaction_parameters += kinetic_law_sbml.getNumParameters()
            num_reaction_parameters += kinetic_law_sbml.getNumLocalParameters()

    model = get_or_create_object(session, Model, id=id)
    model.label = label
    model.name = name
    model.type = type
    model.compartments = sbml_model.getNumCompartments()
    model.species = sbml_model.getNumSpecies()
    model.rules = sbml_model.getNumRules()
    model.reactions = sbml_model.getNumReactions()
    model.global_parameters = sbml_model.getNumParameters()
    model.reaction_parameters = num_reaction_parameters
    model.curated = curated
    model.annotations.extend(annotations)

    session.add(model)

    return model


def parse_model_type(doc, annotations):
    """
    Args:
        doc (:obj:`libsbml.SBMLDocument`): SBML document
        annotations (:obj:`list`: of :obj:`Annotation`): list of annotations

    Returns:
        :obj:`str`: model type
    """

    model = doc.getModel()

    if doc.getPackageRequired('spatial'):
        return 'spatial'

    if doc.getPackageRequired('qual'):
        return 'logical'

    if doc.getPackageRequired('multi'):
        return 'rule-based'

    for annotation in annotations:
        if annotation.namespace == 'mamo' and annotation.id == 'MAMO_0000046':
            return 'ordinary differential equation'

    if doc.getPackageRequired('fbc'):
        return 'flux balance analysis'

    reactions = model.getListOfReactions()
    for i_reaction in range(model.getNumReactions()):
        reaction = reactions.get(i_reaction)
        kinetic_law = reaction.getKineticLaw()
        if kinetic_law:
            has_lower_bound = False
            has_upper_bound = False
            has_flux_value = False
            has_obj_coeff = False

            parameters = kinetic_law.getListOfParameters()
            for i_parameter in range(kinetic_law.getNumParameters()):
                parameter = parameters.get(i_parameter)
                id = parameter.getId()
                if id == 'LOWER_BOUND':
                    has_lower_bound = True
                elif id == 'UPPER_BOUND':
                    has_upper_bound = True
                elif id == 'FLUX_VALUE':
                    has_flux_value = True
                elif id == 'OBJECTIVE_COEFFICIENT':
                    has_obj_coeff = True

            parameters = kinetic_law.getListOfLocalParameters()
            for i_parameter in range(kinetic_law.getNumLocalParameters()):
                parameter = parameters.get(i_parameter)
                id = parameter.getId()
                if id == 'LOWER_BOUND':
                    has_lower_bound = True
                elif id == 'UPPER_BOUND':
                    has_upper_bound = True
                elif id == 'FLUX_VALUE':
                    has_flux_value = True
                elif id == 'OBJECTIVE_COEFFICIENT':
                    has_obj_coeff = True

            if has_lower_bound and has_upper_bound and has_flux_value and has_obj_coeff:
                return 'flux balance analysis'

    return None


def parse_model_annotations(model, session):
    """
    Args:
        model (:obj:`libsbml.Model`): model
        session (:obj:`sqlalchemy.orm.session.Session`): sqlalchemy session

    Returns:
        :obj:`list` of :obj:`Annotation`: list of annotations
    """
    if not model.isSetAnnotation():
        return {}

    annotations_sbml = model.getAnnotation().getChild('RDF').getChild('Description')

    tags = {}
    annotations = []
    attr = libsbml.XMLTriple('resource', 'http://www.w3.org/1999/02/22-rdf-syntax-ns#', 'rdf')
    for i_child in range(annotations_sbml.getNumChildren()):
        child = annotations_sbml.getChild(i_child)

        relationship = child.getName()
        if relationship in ['creator', 'created', 'modified']:
            continue

        for i_bag in range(child.getNumChildren()):
            bag = child.getChild(i_bag)
            if bag.getName() != 'Bag':
                raise ValueError('Expected Bag, got {0}.{1} for model {2}'.format(child.getName(), bag.getName(), model.getId()))

            for i_li in range(bag.getNumChildren()):
                li = bag.getChild(i_li)
                if li.getName() != 'li':
                    raise ValueError('Expected {0}.{1}.li, got {0}.{1}.{2} for model {3}'.format(
                        child.getName(), bag.getName(), li.getName(), model.getId()))

                resource = li.getAttrValue(attr)
                if resource.startswith('http://identifiers.org/'):
                    tmp = resource.split('/')
                    namespace = tmp[3]
                    id = '/'.join(tmp[4:])
                else:
                    namespace = 'url'
                    id = resource
                annotations.append(get_or_create_object(session, Annotation, namespace=namespace, id=id, relationship=relationship))

    return annotations


def get_or_create_object(session, cls, **kwargs):
    """
    Args:
        session (:obj:`sqlalchemy.orm.session.Session`): sqlalchemy session
        cls (:obj:`type`): class to search or create
        **kwargs (:obj:`dict`): dictionary of keyword arguments to pass to filter_by and the class constructor

    Returns:
        :obj:`Base`
    """
    q = session.query(cls).filter_by(**kwargs)
    if q.count():
        return q.first()
    return cls(**kwargs)


def model_to_str(model):
    """
    Args:
        model (:obj:`Model`): model

    Returns:
        :obj:`str`: string representation of model
    """
    str = model.id
    for annotation in model.annotations:
        str += '\n  {}: {}:{}'.format(annotation.relationship, annotation.namespace, annotation.id)
    return str


def export_annotations_to_excel():
    if os.path.isfile(ANNOTATIONS_EXCEL_FILENAME):
        return

    session = get_database_session()

    wb = wc_utils.workbook.core.Workbook()
    ws_models = wb['Models'] = wc_utils.workbook.core.Worksheet()
    ws_model_annotations = wb['Model annotations'] = wc_utils.workbook.core.Worksheet()
    ws_annotations = wb['Annotations'] = wc_utils.workbook.core.Worksheet()
    ws_namespaces = wb['Namespaces'] = wc_utils.workbook.core.Worksheet()

    style = wc_utils.workbook.io.WorkbookStyle()
    style['Models'] = wc_utils.workbook.io.WorksheetStyle(
        head_rows=1, head_columns=1,
        head_row_font_bold=True, head_row_fill_fgcolor='CCCCCC', row_height=15)
    style['Model annotations'] = wc_utils.workbook.io.WorksheetStyle(
        head_rows=1, head_columns=1,
        head_row_font_bold=True, head_row_fill_fgcolor='CCCCCC', row_height=15)
    style['Annotations'] = wc_utils.workbook.io.WorksheetStyle(
        head_rows=1, head_columns=1,
        head_row_font_bold=True, head_row_fill_fgcolor='CCCCCC', row_height=15)
    style['Namespaces'] = wc_utils.workbook.io.WorksheetStyle(
        head_rows=1, head_columns=1,
        head_row_font_bold=True, head_row_fill_fgcolor='CCCCCC', row_height=15)

    ws_models.append(wc_utils.workbook.core.Row([
        'ID', 'Label', 'Name', 'Type',
        'Compartments', 'Species', 'Rules', 'Reactions', 'Global parameters', 'Reaction parameters',
        'Superkingdom', 'Kingdom', 'Phylum', 'Species',
        'Is curated', 'Number annotations']))
    ws_model_annotations.append(wc_utils.workbook.core.Row(['Model', 'Relationship', 'Namespace', 'ID', 'Description']))
    ws_annotations.append(wc_utils.workbook.core.Row(['Relationship', 'Namespace', 'Frequency']))
    ws_namespaces.append(wc_utils.workbook.core.Row(['Namespace', 'Frequency']))

    bio_portal = kinetic_datanator.data_source.bio_portal.BioPortal()
    bio_portal_ontologies = bio_portal.get_ontologies()
    del(bio_portal_ontologies['MAMO'])  # remove MAMO becuse OWL can't be parsed by pronto
    kegg = bioservices.kegg.KEGG()
    reactome = bioservices.reactome.Reactome()
    loaded_ontologies = {}
    ncbi_taxa = ete3.NCBITaxa()
    n_model = session.query(Model).count()
    print('Annotating models ...')
    for i_model, model in enumerate(session.query(Model).order_by(Model.id).all()):
        if i_model % 100 == 0:
            print('  Annotating model {} of {}'.format(i_model + 1, n_model))
        species_name = None
        phylum_name = None
        kingdom_name = None
        superkingdom_name = None
        taxon_id = next((int(float(a.id)) for a in model.annotations if a.namespace == 'taxonomy'), None)
        if taxon_id:
            for taxon_id, rank in ncbi_taxa.get_rank(ncbi_taxa.get_lineage(taxon_id)).items():
                if rank == 'species':
                    species_name = ncbi_taxa.translate_to_names([taxon_id])[0]
                if rank == 'phylum':
                    phylum_name = ncbi_taxa.translate_to_names([taxon_id])[0]
                if rank == 'kingdom':
                    kingdom_name = ncbi_taxa.translate_to_names([taxon_id])[0]
                if rank == 'superkingdom':
                    superkingdom_name = ncbi_taxa.translate_to_names([taxon_id])[0]

        ws_models.append(wc_utils.workbook.core.Row([
            model.id, model.label, model.name, model.type,
            model.compartments or None, model.species or None, model.rules or None, model.reactions or None,
            model.global_parameters or None, model.reaction_parameters or None,
            superkingdom_name, kingdom_name, phylum_name, species_name,
            model.curated, len(model.annotations)
        ]))
        for annotation in sorted(model.annotations, key=lambda ann: (ann.relationship, ann.namespace, ann.id)):

            onto_id = annotation.namespace.upper()
            if onto_id.startswith('OBO.'):
                onto_id = onto_id[4:]
            term_id = annotation.id
            if onto_id in bio_portal_ontologies and term_id.startswith(onto_id + ':'):
                if onto_id not in loaded_ontologies:
                    loaded_ontologies[onto_id] = bio_portal.get_ontology(onto_id)
                onto = loaded_ontologies[onto_id]
                if term_id in onto:
                    description = onto[term_id].name
                else:
                    description = None
            elif annotation.namespace == 'kegg.pathway':
                md = kegg.parse(kegg.get(annotation.id))
                if isinstance(md, dict):
                    description = md['NAME'][0]
                else:
                    description = None
            elif annotation.namespace == 'reactome':
                md = reactome.query_by_id('Pathway', annotation.id)
                if 'displayName' in md:
                    description = md['displayName']
                else:
                    description = None
            elif annotation.namespace == 'taxonomy':
                description = ncbi_taxa.translate_to_names([int(float(annotation.id))])[0]
            else:
                description = None

            ws_model_annotations.append(wc_utils.workbook.core.Row(
                [model.id, annotation.relationship, annotation.namespace, annotation.id, description]))
    print('  done')

    q = session \
        .query(Annotation.relationship, Annotation.namespace, sqlalchemy.func.count(Model._id)) \
        .join(Model, Annotation.models) \
        .group_by(Annotation.relationship, Annotation.namespace) \
        .order_by(sqlalchemy.func.count(Model._id).desc())
    for relationship, namespace, count in q.all():
        ws_annotations.append(wc_utils.workbook.core.Row([relationship, namespace, count]))

    q = session \
        .query(Annotation.namespace, sqlalchemy.func.count(Model._id)) \
        .join(Model, Annotation.models) \
        .group_by(Annotation.namespace) \
        .order_by(sqlalchemy.func.count(Model._id).desc())
    for namespace, count in q.all():
        ws_namespaces.append(wc_utils.workbook.core.Row([namespace, count]))

    wc_utils.workbook.io.ExcelWriter(ANNOTATIONS_EXCEL_FILENAME).run(wb, style=style)


def summarize_models():
    wb = wc_utils.workbook.core.Workbook()
    style = wc_utils.workbook.io.WorkbookStyle()

    ws = wb['Pathways'] = wc_utils.workbook.core.Worksheet()
    style['Pathways'] = wc_utils.workbook.io.WorksheetStyle(
        head_rows=1, head_columns=1, head_row_font_bold=True, head_row_fill_fgcolor='CCCCCC', row_height=15)
    summarize_models_by_pathway(ws)

    ws_species = wb['Species'] = wc_utils.workbook.core.Worksheet()
    ws_phyla = wb['Phyla'] = wc_utils.workbook.core.Worksheet()
    style['Species'] = wc_utils.workbook.io.WorksheetStyle(
        head_rows=1, head_columns=1, head_row_font_bold=True, head_row_fill_fgcolor='CCCCCC', row_height=15)
    style['Phyla'] = wc_utils.workbook.io.WorksheetStyle(
        head_rows=1, head_columns=1, head_row_font_bold=True, head_row_fill_fgcolor='CCCCCC', row_height=15)
    summarize_models_by_taxonomy(ws_species, ws_phyla)

    ws = wb['Mathematical types'] = wc_utils.workbook.core.Worksheet()
    style['Mathematical types'] = wc_utils.workbook.io.WorksheetStyle(
        head_rows=1, head_columns=1, head_row_font_bold=True, head_row_fill_fgcolor='CCCCCC', row_height=15)
    summarize_models_by_mathematical_type(ws)

    wc_utils.workbook.io.ExcelWriter(SUMMARY_EXCEL_FILENAME).run(wb, style=style)


def summarize_models_by_pathway(ws):
    """
    Args:
        wc_utils.workbook.core.Worksheet
    """
    session = get_database_session()
    #bio_portal = kinetic_datanator.data_source.bio_portal.BioPortal()
    #onto = bio_portal.get_ontology('EFO')
    #self.assertEqual(onto['SBO:0000001'].name, 'rate law')


def summarize_models_by_taxonomy(ws_species, ws_phyla):
    """
    Args:
        wc_utils.workbook.core.Worksheet
    """
    session = get_database_session()
    q_annotated = session \
        .query(Annotation.id, Model.curated, sqlalchemy.func.count(Model._id)) \
        .join(Model, Annotation.models) \
        .filter(Annotation.namespace == 'taxonomy') \
        .group_by(Annotation.id, Model.curated)

    annotated_model_ids = [m[0] for m in session
                           .query(Model._id)
                           .join(Annotation, Model.annotations)
                           .filter(Annotation.namespace == 'taxonomy')
                           .group_by(Model._id)
                           .all()]
    q_unannotated = session \
        .query(Model.curated, sqlalchemy.func.count(Model._id)) \
        .filter(~Model._id.in_(annotated_model_ids)) \
        .group_by(Model.curated)
    count_unannotated = {}
    for curated, count in q_unannotated.all():
        count_unannotated[curated] = count

    ncbi_taxa = ete3.NCBITaxa()

    species = {}
    phyla = {}
    for model_taxon_id, model_curated, count in q_annotated.all():
        model_taxon_id = int(float(model_taxon_id))
        species_name = None
        phylum_name = None
        kingdom_name = None
        superkingdom_name = None
        for taxon_id, rank in ncbi_taxa.get_rank(ncbi_taxa.get_lineage(model_taxon_id)).items():
            if rank == 'species':
                species_name = ncbi_taxa.translate_to_names([taxon_id])[0]
            if rank == 'phylum':
                phylum_name = ncbi_taxa.translate_to_names([taxon_id])[0]
            if rank == 'kingdom':
                kingdom_name = ncbi_taxa.translate_to_names([taxon_id])[0]
            if rank == 'superkingdom':
                superkingdom_name = ncbi_taxa.translate_to_names([taxon_id])[0]

        if (superkingdom_name, kingdom_name, phylum_name, species_name) not in species:
            species[(superkingdom_name, kingdom_name, phylum_name, species_name)] = {True: 0, False: 0}
        species[(superkingdom_name, kingdom_name, phylum_name, species_name)][model_curated] += count

        if (superkingdom_name, kingdom_name, phylum_name) not in phyla:
            phyla[(superkingdom_name, kingdom_name, phylum_name)] = {True: 0, False: 0}
        phyla[(superkingdom_name, kingdom_name, phylum_name)][model_curated] += count

    for (superkingdom_name, kingdom_name, phylum_name, species_name), counts in species.items():
        ws_species.append(wc_utils.workbook.core.Row([
            species_name or '<Annotated rank above species>',
            phylum_name,
            kingdom_name,
            superkingdom_name,
            counts[True] or None,
            counts[False] or None]))
    ws_species.sort(key=lambda row: (row[-2] or 0) + (row[-1] or 0), reverse=True)
    ws_species.insert(0, wc_utils.workbook.core.Row(['Not annotated', None, None, count_unannotated[True], count_unannotated[False]]))
    ws_species.insert(0, wc_utils.workbook.core.Row(['Species', 'Phylum', 'Kingdom', 'Superkingdom', 'Curated', 'Non-curated']))

    for (superkingdom_name, kingdom_name, phylum_name), counts in phyla.items():
        ws_phyla.append(wc_utils.workbook.core.Row([
            superkingdom_name,
            kingdom_name,
            phylum_name or '<Annotated rank above phylum>',
            counts[True] or None,
            counts[False] or None]))
    ws_phyla.sort(key=lambda row: (row[-2] or 0) + (row[-1] or 0), reverse=True)
    ws_phyla.insert(0, wc_utils.workbook.core.Row(['Not annotated', None, count_unannotated[True], count_unannotated[False]]))
    ws_phyla.insert(0, wc_utils.workbook.core.Row(['Superkingdom', 'Kingdom', 'Phylum', 'Curated', 'Non-curated']))


def summarize_models_by_mathematical_type(ws):
    """
    Args:
        wc_utils.workbook.core.Worksheet
    """
    session = get_database_session()

    q = session.query(Model.type, Model.curated, sqlalchemy.func.count(Model._id)) \
        .group_by(Model.type, Model.curated) \
        .order_by(Model.type)
    data = {}
    for type, curated, count in q.all():
        if type not in data:
            data[type] = {}
        data[type][curated] = count

    ws.append(wc_utils.workbook.core.Row(['Type', 'Curated', 'Non-curated']))
    for type in data.keys():
        ws.append(wc_utils.workbook.core.Row([
            type[0].upper() + type[1:] if type else 'Unknown',
            data[type][True] if True in data[type] else None,
            data[type][False] if False in data[type] else None,
        ]))

Base = sqlalchemy.ext.declarative.declarative_base()
# :obj:`Base`: base model for local sqlite database

model_annotation = sqlalchemy.Table(
    'model_annotation', Base.metadata,
    sqlalchemy.Column('model__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('model._id'), index=True),
    sqlalchemy.Column('annotation__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('annotation._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Model:Annotation many-to-many association table


class Model(Base):
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    id = sqlalchemy.Column(sqlalchemy.String())
    label = sqlalchemy.Column(sqlalchemy.String())
    name = sqlalchemy.Column(sqlalchemy.String())
    type = sqlalchemy.Column(sqlalchemy.String())
    curated = sqlalchemy.Column(sqlalchemy.Boolean())
    compartments = sqlalchemy.Column(sqlalchemy.Integer())
    species = sqlalchemy.Column(sqlalchemy.Integer())
    rules = sqlalchemy.Column(sqlalchemy.Integer())
    reactions = sqlalchemy.Column(sqlalchemy.Integer())
    global_parameters = sqlalchemy.Column(sqlalchemy.Integer())
    reaction_parameters = sqlalchemy.Column(sqlalchemy.Integer())
    annotations = sqlalchemy.orm.relationship('Annotation', secondary=model_annotation, backref=sqlalchemy.orm.backref('models'))

    __tablename__ = 'model'


class Annotation(Base):
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    namespace = sqlalchemy.Column(sqlalchemy.String(), index=True)
    id = sqlalchemy.Column(sqlalchemy.String(), index=True)
    relationship = sqlalchemy.Column(sqlalchemy.String(), index=True)

    __tablename__ = 'annotation'


if __name__ == "__main__":
    create_data_directory()
    download_biomodels()
    setup_database()
    load_database()
    export_annotations_to_excel()
    summarize_models()
