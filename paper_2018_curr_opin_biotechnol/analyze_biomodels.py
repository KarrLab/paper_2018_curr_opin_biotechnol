""" Analysis of the content of BioModels

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-31
:Copyright: 2017, Karr Lab
:License: MIT
"""

import glob
import libsbml
import os
import sqlalchemy
import sqlalchemy.ext.declarative
import sqlalchemy.orm
import subprocess
import wc_utils.workbook.core
import wc_utils.workbook.io

RELEASE_NAME = 'BioModels_Database-r30_pub-sbml_files'
SBML_FILES_URL = 'ftp://ftp.ebi.ac.uk/pub/databases/biomodels/releases/2016-05-10/{}.tar.bz2'.format(RELEASE_NAME)
DATA_DIRNAME = os.path.join(os.path.dirname(__file__), 'data')
SBML_FILES_ARCHIVE_FILENAME = os.path.join(DATA_DIRNAME, 'sbml_files.tar.bz2')
SBML_FILES_DIRNAME = os.path.join(DATA_DIRNAME, 'sbml_files')
SBML_FILES_DATABASE_FILENAME = os.path.join(DATA_DIRNAME, 'sbml_files.sqlite')
ANNOTATIONS_EXCEL_FILENAME = os.path.join(DATA_DIRNAME, 'models.xlsx')


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
    model.type = parse_model_type(doc)
    model.compartments = sbml_model.getNumCompartments()
    model.species = sbml_model.getNumSpecies()
    model.rules = sbml_model.getNumRules()
    model.reactions = sbml_model.getNumReactions()
    model.global_parameters = sbml_model.getNumParameters()
    model.reaction_parameters = num_reaction_parameters
    model.curated = curated
    model.annotations.extend(parse_model_annotations(sbml_model, session))

    session.add(model)

    return model


def parse_model_type(doc):
    """
    Args:
        doc (:obj:`libsbml.SBMLDocument`): SBML document

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
    session = get_database_session()

    wb = wc_utils.workbook.core.Workbook()
    ws_mod = wb['Models'] = wc_utils.workbook.core.Worksheet()
    ws_ann = wb['Annotations'] = wc_utils.workbook.core.Worksheet()

    ws_mod.append(wc_utils.workbook.core.Row([
        'ID', 'Label', 'Name', 'Type',
        'Compartments', 'Species', 'Rules', 'Reactions', 'Global parameters', 'Reaction parameters',
        'Is curated', 'Number annotations']))
    ws_ann.append(wc_utils.workbook.core.Row(['Model', 'Relationship', 'Namespace', 'ID']))

    for model in session.query(Model).all():
        ws_mod.append(wc_utils.workbook.core.Row([
            model.id, model.label, model.name, model.type,
            model.compartments, model.species, model.rules, model.reactions, model.global_parameters, model.reaction_parameters,
            model.curated, len(model.annotations)
        ]))
        for annotation in model.annotations:
            ws_ann.append(wc_utils.workbook.core.Row([model.id, annotation.relationship, annotation.namespace, annotation.id]))

    wc_utils.workbook.io.ExcelWriter(ANNOTATIONS_EXCEL_FILENAME).run(wb)

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
