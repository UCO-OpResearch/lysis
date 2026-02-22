"""Comprehensive tests for the lysis.data.dataspec module.

Tests cover:
- parse_shape() with integers, variable dims, string references, and errors
- check_dataset_spec() with dtype casting, fixed/variable dims, and params
- DataSetSpec and DataCollectionSpec frozen dataclass behavior
- The dataspec registry dict (versions, tags, collections, dataset names)
- The tags dict (aliases and resolution)
"""

import copy
from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from lysis.config.constants import CONST
from lysis.data.dataspec import (
    DataCollectionSpec,
    DataSetSpec,
    DataSpec,
    check_dataset_spec,
    dataspec,
    parse_shape,
    tags,
)


# ---------------------------------------------------------------------------
# parse_shape tests
# ---------------------------------------------------------------------------
class TestParseShape:
    """Tests for the parse_shape function."""

    def test_integers_only(self):
        """Integer-only shape is returned unchanged."""
        assert parse_shape((100, 3)) == (100, 3)

    def test_variable_dimension(self):
        """-1 sentinel passes through unchanged."""
        assert parse_shape((-1, 5)) == (-1, 5)

    def test_string_reference_resolves(self, sample_params):
        """String element is resolved via the params dict."""
        assert parse_shape((-1, "macro_params.total_molecules"), sample_params) == (
            -1,
            10,
        )

    def test_mixed_int_and_string(self, sample_params):
        """Mixed integers and strings are handled together."""
        assert parse_shape((3, "macro_params.rows"), sample_params) == (3, 5)

    def test_invalid_type_raises_runtime_error(self):
        """Non-int, non-str element raises RuntimeError."""
        with pytest.raises(RuntimeError, match="Incorrect shape format"):
            parse_shape((3.14,))

    def test_string_without_params_raises(self):
        """String element with params=None raises TypeError or KeyError."""
        with pytest.raises((TypeError, KeyError)):
            parse_shape(("macro_params.x",), None)

    def test_empty_shape(self):
        """Empty shape tuple returns empty tuple."""
        assert parse_shape(()) == ()

    def test_single_integer(self):
        """Single-element integer shape is handled."""
        assert parse_shape((42,)) == (42,)


# ---------------------------------------------------------------------------
# check_dataset_spec tests
# ---------------------------------------------------------------------------
class TestCheckDatasetSpec:
    """Tests for the check_dataset_spec function."""

    @pytest.mark.parametrize(
        "data_dtype, spec_dtype, expected",
        [
            pytest.param(np.float64, np.float64, True, id="same_dtype"),
            pytest.param(np.int32, np.float64, True, id="int32_to_float64_castable"),
            pytest.param(
                np.float64, np.int32, False, id="float64_to_int32_not_castable"
            ),
        ],
    )
    def test_dtype_casting(self, data_dtype, spec_dtype, expected):
        """Dtype compatibility is checked via np.can_cast."""
        data = np.zeros((5,), dtype=data_dtype)
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
            dtype=spec_dtype,
            shape=(-1,),
        )
        assert check_dataset_spec(data, spec) is expected

    def test_fixed_shape_match(self):
        """Fixed shape that exactly matches returns True."""
        data = np.zeros((3, 4), dtype=np.float64)
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
            dtype=np.float64,
            shape=(3, 4),
        )
        assert check_dataset_spec(data, spec) is True

    def test_fixed_shape_mismatch(self):
        """Fixed shape that does not match returns False."""
        data = np.zeros((3, 5), dtype=np.float64)
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
            dtype=np.float64,
            shape=(3, 4),
        )
        assert check_dataset_spec(data, spec) is False

    def test_variable_dim_accepts_any_size(self):
        """-1 dimension accepts any actual size."""
        data = np.zeros((99,), dtype=np.float64)
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
            dtype=np.float64,
            shape=(-1,),
        )
        assert check_dataset_spec(data, spec) is True

    def test_dynamic_shape_from_params(self, sample_params):
        """String shape resolved from params matches the data."""
        data = np.zeros((10,), dtype=np.float64)
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
            dtype=np.float64,
            shape=("macro_params.total_molecules",),
        )
        assert check_dataset_spec(data, spec, params=sample_params) is True

    def test_mixed_variable_and_fixed_match(self):
        """-1 with a fixed dim: data matches when fixed dim is correct."""
        data = np.zeros((7, 10), dtype=np.float64)
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
            dtype=np.float64,
            shape=(-1, 10),
        )
        assert check_dataset_spec(data, spec) is True

    def test_mixed_variable_and_fixed_mismatch(self):
        """-1 with a fixed dim: data fails when fixed dim is wrong."""
        data = np.zeros((7, 9), dtype=np.float64)
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
            dtype=np.float64,
            shape=(-1, 10),
        )
        assert check_dataset_spec(data, spec) is False


# ---------------------------------------------------------------------------
# DataSetSpec frozen dataclass tests
# ---------------------------------------------------------------------------
class TestDataSetSpec:
    """Tests for the DataSetSpec frozen dataclass."""

    def test_frozen_immutability(self):
        """Assigning to a field on a frozen DataSetSpec raises FrozenInstanceError."""
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
            dtype=np.float64,
        )
        with pytest.raises(FrozenInstanceError):
            spec.dtype = np.int32

    def test_default_shape(self):
        """Default shape is (-1,)."""
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
            dtype=np.float64,
        )
        assert spec.shape == (-1,)

    @pytest.mark.parametrize(
        "storage_type",
        [
            CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
            CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
            CONST.DATASET_STORAGE_TYPE.FILE_JSON,
            CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
            CONST.DATASET_STORAGE_TYPE.HDF5_ATTR,
        ],
        ids=["FILE_TEXT", "FILE_BINARY", "FILE_JSON", "HDF5_DATASET", "HDF5_ATTR"],
    )
    def test_all_storage_types_constructable(self, storage_type):
        """Each of the five storage types can construct a DataSetSpec."""
        spec = DataSetSpec(
            dataset_storage_type=storage_type,
            dtype=np.float64,
        )
        assert spec.dataset_storage_type is storage_type


# ---------------------------------------------------------------------------
# DataCollectionSpec frozen dataclass tests
# ---------------------------------------------------------------------------
class TestDataCollectionSpec:
    """Tests for the DataCollectionSpec frozen dataclass."""

    def test_frozen_immutability(self):
        """Assigning to a field on a frozen DataCollectionSpec raises FrozenInstanceError."""
        spec = DataCollectionSpec(
            simulations_combined=True,
            params=None,
            data={},
        )
        with pytest.raises(FrozenInstanceError):
            spec.simulations_combined = False

    def test_params_none(self):
        """DataCollectionSpec can be created with params=None."""
        spec = DataCollectionSpec(
            simulations_combined=True,
            params=None,
            data={},
        )
        assert spec.params is None

    @pytest.mark.parametrize(
        "combined", [True, False], ids=["combined", "per_sim"]
    )
    def test_simulations_combined_flag(self, combined):
        """simulations_combined can be True or False."""
        spec = DataCollectionSpec(
            simulations_combined=combined,
            params=None,
            data={},
        )
        assert spec.simulations_combined is combined


# ---------------------------------------------------------------------------
# dataspec registry tests
# ---------------------------------------------------------------------------
class TestDataspecRegistry:
    """Tests for the dataspec module-level registry dict."""

    @pytest.mark.parametrize("version", ["v1.99.0", "v2.0.0"])
    def test_version_keys_exist(self, version):
        """All canonical version keys are present."""
        assert version in dataspec

    @pytest.mark.parametrize(
        "tag, version",
        [
            ("fortran", "v1.99.0"),
            ("hdf5", "v2.0.0"),
        ],
    )
    def test_tag_aliases_resolve(self, tag, version):
        """Tag aliases point to the same dict object as the canonical version."""
        assert dataspec[tag] is dataspec[version]

    def test_current_tag_resolves_to_hdf5(self):
        """The 'current' tag resolves to the same data as 'hdf5'."""
        assert dataspec["current"] is dataspec["hdf5"]

    @pytest.mark.parametrize("version", ["v1.99.0", "v2.0.0"])
    @pytest.mark.parametrize(
        "collection", ["microscale_out", "macroscale_in", "macroscale_out"]
    )
    def test_each_version_has_expected_collections(self, version, collection):
        """Each version has all three expected collection names."""
        assert collection in dataspec[version]

    @pytest.mark.parametrize("version", ["v1.99.0", "v2.0.0"])
    def test_collection_datasets_have_valid_dtype(self, version):
        """Every dataset in every collection has a non-None dtype."""
        for coll_name, coll_spec in dataspec[version].items():
            for ds_name, ds_spec in coll_spec.data.items():
                assert ds_spec.dtype is not None, (
                    f"{version}/{coll_name}/{ds_name} has dtype=None"
                )

    def test_v1_microscale_out_dataset_names(self):
        """v1.99.0 microscale_out has the expected legacy dataset names."""
        expected = {
            "micro_log",
            "firstPLi",
            "lasttPA",
            "lyscomplete",
            "lysis",
            "PLi",
            "tPA_time",
            "tPAPLiunbd",
            "tPAunbind",
        }
        actual = set(dataspec["v1.99.0"]["microscale_out"].data.keys())
        assert actual == expected

    def test_v2_microscale_out_dataset_names(self):
        """v2.0.0 microscale_out has the expected modern dataset names."""
        expected = {
            "micro_log",
            "pli_first_time",
            "tpa_final_num",
            "fiber_degraded",
            "sim_final_time",
            "pli_generated_num",
            "tpa_leaving_time",
            "tpa_unbound_by_pli",
            "tpa_unbound_kinetic",
        }
        actual = set(dataspec["v2.0.0"]["microscale_out"].data.keys())
        assert actual == expected


# ---------------------------------------------------------------------------
# tags dict tests
# ---------------------------------------------------------------------------
class TestTags:
    """Tests for the tags module-level dict."""

    @pytest.mark.parametrize("tag", ["fortran", "hdf5", "current"])
    def test_expected_tags_exist(self, tag):
        """All expected tag names are present."""
        assert tag in tags

    def test_fortran_maps_to_v1(self):
        """'fortran' tag maps to 'v1.99.0'."""
        assert tags["fortran"] == "v1.99.0"

    def test_hdf5_maps_to_v2(self):
        """'hdf5' tag maps to 'v2.0.0'."""
        assert tags["hdf5"] == "v2.0.0"

    def test_current_maps_to_hdf5(self):
        """'current' tag maps to 'hdf5'."""
        assert tags["current"] == "hdf5"


# ---------------------------------------------------------------------------
# DataSpec wrapper class tests
# ---------------------------------------------------------------------------
class TestDataSpec:
    """Tests for the DataSpec wrapper class."""

    def test_version_property(self):
        """DataSpec.version returns the version string."""
        spec = dataspec["v2.0.0"]
        assert spec.version == "v2.0.0"

    def test_v199_version_property(self):
        """DataSpec.version returns the version string for v1.99.0."""
        spec = dataspec["v1.99.0"]
        assert spec.version == "v1.99.0"

    def test_getitem(self):
        """DataSpec supports dict-style key access."""
        spec = dataspec["v2.0.0"]
        coll = spec["microscale_out"]
        assert isinstance(coll, DataCollectionSpec)

    def test_missing_key_raises(self):
        """Accessing a nonexistent key raises KeyError."""
        spec = dataspec["v2.0.0"]
        with pytest.raises(KeyError):
            spec["nonexistent"]

    def test_contains(self):
        """DataSpec supports 'in' operator."""
        spec = dataspec["v2.0.0"]
        assert "microscale_out" in spec
        assert "nonexistent" not in spec

    def test_iter(self):
        """DataSpec is iterable over collection names."""
        spec = dataspec["v2.0.0"]
        names = list(spec)
        assert "microscale_out" in names
        assert "macroscale_out" in names

    def test_len(self):
        """DataSpec has a length equal to the number of collections."""
        spec = dataspec["v2.0.0"]
        assert len(spec) == 3

    def test_items(self):
        """DataSpec.items() yields (name, DataCollectionSpec) pairs."""
        spec = dataspec["v2.0.0"]
        items = dict(spec.items())
        assert "microscale_out" in items
        assert isinstance(items["microscale_out"], DataCollectionSpec)

    def test_keys(self):
        """DataSpec.keys() returns collection names."""
        spec = dataspec["v2.0.0"]
        assert set(spec.keys()) == {"microscale_out", "macroscale_in", "macroscale_out"}

    def test_values(self):
        """DataSpec.values() returns DataCollectionSpec objects."""
        spec = dataspec["v2.0.0"]
        for v in spec.values():
            assert isinstance(v, DataCollectionSpec)

    def test_repr(self):
        """DataSpec repr includes version and collection names."""
        spec = dataspec["v2.0.0"]
        r = repr(spec)
        assert "v2.0.0" in r
        assert "microscale_out" in r

    def test_tag_alias_shares_object(self):
        """Tag aliases point to the same DataSpec object."""
        assert dataspec["hdf5"] is dataspec["v2.0.0"]
        assert dataspec["fortran"] is dataspec["v1.99.0"]

    def test_tag_alias_version_is_canonical(self):
        """Tag-aliased DataSpec has the canonical version, not the tag name."""
        assert dataspec["hdf5"].version == "v2.0.0"
        assert dataspec["current"].version == "v2.0.0"


# ---------------------------------------------------------------------------
# Hidden field population tests
# ---------------------------------------------------------------------------
class TestHiddenFields:
    """Tests for auto-populated hidden fields on DataCollectionSpec and DataSetSpec."""

    def test_collection_spec_has_version(self):
        """DataCollectionSpec.version is set by DataSpec."""
        coll = dataspec["v2.0.0"]["microscale_out"]
        assert coll.version == "v2.0.0"

    def test_collection_spec_has_collection(self):
        """DataCollectionSpec.collection is set to its own name."""
        coll = dataspec["v2.0.0"]["microscale_out"]
        assert coll.collection == "microscale_out"

    def test_dataset_spec_has_version(self):
        """DataSetSpec.version is set by DataSpec."""
        ds = dataspec["v2.0.0"]["microscale_out"].data["pli_first_time"]
        assert ds.version == "v2.0.0"

    def test_dataset_spec_has_collection(self):
        """DataSetSpec.collection is set to the parent collection name."""
        ds = dataspec["v2.0.0"]["microscale_out"].data["pli_first_time"]
        assert ds.collection == "microscale_out"

    def test_dataset_spec_has_name(self):
        """DataSetSpec.name is set to the dataset key."""
        ds = dataspec["v2.0.0"]["microscale_out"].data["pli_first_time"]
        assert ds.name == "pli_first_time"

    def test_params_spec_has_fields(self):
        """Params DataSetSpec has version, collection, and name='params'."""
        params = dataspec["v2.0.0"]["microscale_out"].params
        assert params.version == "v2.0.0"
        assert params.collection == "microscale_out"
        assert params.name == "params"

    def test_standalone_spec_defaults_empty(self):
        """DataSetSpec created outside DataSpec has empty hidden fields."""
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
            dtype=np.float64,
        )
        assert spec.version == ""
        assert spec.collection == ""
        assert spec.name == ""

    def test_standalone_collection_spec_defaults_empty(self):
        """DataCollectionSpec created outside DataSpec has empty hidden fields."""
        coll = DataCollectionSpec(
            simulations_combined=True,
            params=None,
            data={},
        )
        assert coll.version == ""
        assert coll.collection == ""

    def test_fields_not_in_constructor(self):
        """Hidden fields cannot be passed as constructor arguments."""
        with pytest.raises(TypeError):
            DataSetSpec(
                dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                dtype=np.float64,
                name="should_fail",
            )

    def test_v199_fields_populated(self):
        """Hidden fields work for v1.99.0 specs too."""
        coll = dataspec["v1.99.0"]["macroscale_out"]
        assert coll.version == "v1.99.0"
        assert coll.collection == "macroscale_out"
        ds = coll.data["tsave"]
        assert ds.version == "v1.99.0"
        assert ds.collection == "macroscale_out"
        assert ds.name == "tsave"


# ---------------------------------------------------------------------------
# DataSetSpec.__copy__ and __deepcopy__ tests
# ---------------------------------------------------------------------------
class TestDataSetSpecCopy:
    """Tests for DataSetSpec.__copy__ and __deepcopy__."""

    def test_copy_returns_different_object(self, text_spec):
        """copy.copy creates a distinct object."""
        copied = copy.copy(text_spec)
        assert copied is not text_spec

    def test_copy_preserves_init_fields(self, binary_spec):
        """copy.copy preserves all init field values."""
        copied = copy.copy(binary_spec)
        assert copied.dataset_storage_type == binary_spec.dataset_storage_type
        assert copied.dtype == binary_spec.dtype
        assert copied.data_location == binary_spec.data_location
        assert copied.shape == binary_spec.shape
        assert copied.delimiter == binary_spec.delimiter

    def test_copy_preserves_empty_hidden_fields(self, text_spec):
        """copy.copy preserves default (empty) hidden fields."""
        copied = copy.copy(text_spec)
        assert copied.version == ""
        assert copied.collection == ""
        assert copied.name == ""

    def test_copy_preserves_populated_hidden_fields(self):
        """copy.copy preserves hidden fields populated by the dataspec registry."""
        spec = dataspec["v2.0.0"]["microscale_out"].data["pli_first_time"]
        copied = copy.copy(spec)
        assert copied.version == "v2.0.0"
        assert copied.collection == "microscale_out"
        assert copied.name == "pli_first_time"

    def test_copy_result_is_frozen(self, text_spec):
        """Copied DataSetSpec is still frozen."""
        copied = copy.copy(text_spec)
        with pytest.raises(FrozenInstanceError):
            copied.dtype = np.int32

    def test_deepcopy_returns_different_object(self, text_spec):
        """copy.deepcopy creates a distinct object."""
        deep = copy.deepcopy(text_spec)
        assert deep is not text_spec

    def test_deepcopy_preserves_init_fields(self, binary_spec):
        """copy.deepcopy preserves all init field values."""
        deep = copy.deepcopy(binary_spec)
        assert deep.dataset_storage_type == binary_spec.dataset_storage_type
        assert deep.dtype == binary_spec.dtype
        assert deep.data_location == binary_spec.data_location
        assert deep.shape == binary_spec.shape
        assert deep.delimiter == binary_spec.delimiter

    def test_deepcopy_preserves_empty_hidden_fields(self, text_spec):
        """copy.deepcopy preserves default (empty) hidden fields."""
        deep = copy.deepcopy(text_spec)
        assert deep.version == ""
        assert deep.collection == ""
        assert deep.name == ""

    def test_deepcopy_preserves_populated_hidden_fields(self):
        """copy.deepcopy preserves hidden fields populated by the dataspec registry."""
        spec = dataspec["v2.0.0"]["microscale_out"].data["pli_first_time"]
        deep = copy.deepcopy(spec)
        assert deep.version == "v2.0.0"
        assert deep.collection == "microscale_out"
        assert deep.name == "pli_first_time"

    def test_deepcopy_result_is_frozen(self, text_spec):
        """Deep-copied DataSetSpec is still frozen."""
        deep = copy.deepcopy(text_spec)
        with pytest.raises(FrozenInstanceError):
            deep.dtype = np.int32

    def test_deepcopy_structured_dtype(self):
        """copy.deepcopy preserves a structured numpy dtype."""
        structured_spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
            dtype=np.dtype([("time", np.float64), ("index", np.int32)]),
            delimiter=",",
        )
        deep = copy.deepcopy(structured_spec)
        assert deep is not structured_spec
        assert deep.dtype == structured_spec.dtype

    def test_deepcopy_memo_registers_object(self, text_spec):
        """copy.deepcopy registers the result in memo under the original object's id."""
        memo = {}
        deep = copy.deepcopy(text_spec, memo)
        assert id(text_spec) in memo
        assert memo[id(text_spec)] is deep


# ---------------------------------------------------------------------------
# DataCollectionSpec.__copy__ and __deepcopy__ tests
# ---------------------------------------------------------------------------
class TestDataCollectionSpecCopy:
    """Tests for DataCollectionSpec.__copy__ and __deepcopy__."""

    def test_copy_returns_different_object(self, combined_collection_spec):
        """copy.copy creates a distinct object."""
        copied = copy.copy(combined_collection_spec)
        assert copied is not combined_collection_spec

    def test_copy_preserves_simulations_combined(self, combined_collection_spec):
        """copy.copy preserves the simulations_combined flag."""
        copied = copy.copy(combined_collection_spec)
        assert copied.simulations_combined == combined_collection_spec.simulations_combined

    def test_copy_shares_params(self, combined_collection_spec):
        """copy.copy shares the params object (shallow copy semantics)."""
        copied = copy.copy(combined_collection_spec)
        assert copied.params is combined_collection_spec.params

    def test_copy_shares_data_dict(self, combined_collection_spec):
        """copy.copy shares the data dict (shallow copy semantics)."""
        copied = copy.copy(combined_collection_spec)
        assert copied.data is combined_collection_spec.data

    def test_copy_preserves_empty_hidden_fields(self, combined_collection_spec):
        """copy.copy preserves default (empty) hidden fields."""
        copied = copy.copy(combined_collection_spec)
        assert copied.version == ""
        assert copied.collection == ""

    def test_copy_preserves_populated_hidden_fields(self):
        """copy.copy preserves hidden fields populated by the dataspec registry."""
        coll = dataspec["v2.0.0"]["microscale_out"]
        copied = copy.copy(coll)
        assert copied.version == "v2.0.0"
        assert copied.collection == "microscale_out"

    def test_copy_with_none_params(self):
        """copy.copy works when params is None."""
        coll = DataCollectionSpec(simulations_combined=True, params=None, data={})
        copied = copy.copy(coll)
        assert copied.params is None
        assert copied.simulations_combined is True

    def test_copy_result_is_frozen(self, combined_collection_spec):
        """Copied DataCollectionSpec is still frozen."""
        copied = copy.copy(combined_collection_spec)
        with pytest.raises(FrozenInstanceError):
            copied.simulations_combined = False

    def test_deepcopy_returns_different_object(self, combined_collection_spec):
        """copy.deepcopy creates a distinct object."""
        deep = copy.deepcopy(combined_collection_spec)
        assert deep is not combined_collection_spec

    def test_deepcopy_creates_new_data_dict(self, combined_collection_spec):
        """copy.deepcopy creates a new data dict object."""
        deep = copy.deepcopy(combined_collection_spec)
        assert deep.data is not combined_collection_spec.data

    def test_deepcopy_creates_new_dataset_specs(self, combined_collection_spec):
        """copy.deepcopy creates new DataSetSpec objects for each entry in data."""
        deep = copy.deepcopy(combined_collection_spec)
        for key in combined_collection_spec.data:
            assert deep.data[key] is not combined_collection_spec.data[key]

    def test_deepcopy_creates_new_params(self, combined_collection_spec):
        """copy.deepcopy creates a new params DataSetSpec object."""
        deep = copy.deepcopy(combined_collection_spec)
        assert deep.params is not combined_collection_spec.params

    def test_deepcopy_preserves_data_values(self, combined_collection_spec):
        """copy.deepcopy preserves the field values of each copied DataSetSpec."""
        deep = copy.deepcopy(combined_collection_spec)
        for key, spec in combined_collection_spec.data.items():
            assert deep.data[key].dtype == spec.dtype
            assert deep.data[key].dataset_storage_type == spec.dataset_storage_type
            assert deep.data[key].data_location == spec.data_location

    def test_deepcopy_preserves_empty_hidden_fields(self, combined_collection_spec):
        """copy.deepcopy preserves default (empty) hidden fields."""
        deep = copy.deepcopy(combined_collection_spec)
        assert deep.version == ""
        assert deep.collection == ""

    def test_deepcopy_preserves_populated_hidden_fields(self):
        """copy.deepcopy preserves hidden fields populated by the dataspec registry."""
        coll = dataspec["v2.0.0"]["microscale_out"]
        deep = copy.deepcopy(coll)
        assert deep.version == "v2.0.0"
        assert deep.collection == "microscale_out"

    def test_deepcopy_propagates_child_hidden_fields(self):
        """DataSetSpec children deep-copied from a populated collection keep their hidden fields."""
        coll = dataspec["v2.0.0"]["microscale_out"]
        deep = copy.deepcopy(coll)
        ds = deep.data["pli_first_time"]
        assert ds.version == "v2.0.0"
        assert ds.collection == "microscale_out"
        assert ds.name == "pli_first_time"

    def test_deepcopy_with_none_params(self):
        """copy.deepcopy works when params is None."""
        coll = DataCollectionSpec(simulations_combined=True, params=None, data={})
        deep = copy.deepcopy(coll)
        assert deep.params is None

    def test_deepcopy_result_is_frozen(self, combined_collection_spec):
        """Deep-copied DataCollectionSpec is still frozen."""
        deep = copy.deepcopy(combined_collection_spec)
        with pytest.raises(FrozenInstanceError):
            deep.simulations_combined = False

    def test_deepcopy_memo_registers_object(self, combined_collection_spec):
        """copy.deepcopy registers the result in memo under the original object's id."""
        memo = {}
        deep = copy.deepcopy(combined_collection_spec, memo)
        assert id(combined_collection_spec) in memo
        assert memo[id(combined_collection_spec)] is deep


# ---------------------------------------------------------------------------
# DataCollectionSpec.replace() tests
# ---------------------------------------------------------------------------
class TestDataCollectionSpecReplace:
    """Tests for DataCollectionSpec.replace().

    The fixture ``combined_collection_spec`` has:
        simulations_combined=True, params=json_spec, data={"arr": text_spec}

    where json_spec has dtype=dict and text_spec has dtype=np.float64.
    """

    # -----------------------------------------------------------------------
    # Top-level field replacements
    # -----------------------------------------------------------------------

    def test_replace_bool_field_directly(self, combined_collection_spec):
        """Replacing simulations_combined with a bool sets the new value."""
        new = combined_collection_spec.replace(simulations_combined=False)
        assert new.simulations_combined is False

    def test_replace_params_with_full_spec(self, combined_collection_spec, hdf5_attr_spec):
        """Replacing params with a DataSetSpec uses that object directly."""
        new = combined_collection_spec.replace(params=hdf5_attr_spec)
        assert new.params is hdf5_attr_spec

    def test_replace_params_with_dict_updates_specified_fields(
        self, combined_collection_spec
    ):
        """Replacing params with a dict updates only the specified fields."""
        new = combined_collection_spec.replace(params={"dtype": np.int32})
        assert new.params.dtype == np.int32
        # Fields not mentioned are preserved from the original params spec.
        assert (
            new.params.dataset_storage_type
            == combined_collection_spec.params.dataset_storage_type
        )
        assert new.params.data_location == combined_collection_spec.params.data_location

    def test_replace_params_with_dict_returns_new_params_object(
        self, combined_collection_spec
    ):
        """Dict-form params replacement creates a new DataSetSpec instance."""
        new = combined_collection_spec.replace(params={"dtype": np.int32})
        assert new.params is not combined_collection_spec.params

    # -----------------------------------------------------------------------
    # Named data-entry replacements
    # -----------------------------------------------------------------------

    def test_replace_data_entry_with_full_spec(
        self, combined_collection_spec, binary_spec
    ):
        """Replacing a named data entry with a DataSetSpec uses it directly."""
        new = combined_collection_spec.replace(arr=binary_spec)
        assert new.data["arr"] is binary_spec

    def test_replace_data_entry_with_dict_updates_specified_fields(
        self, combined_collection_spec
    ):
        """Dict-form data-entry replacement updates only the specified fields."""
        original_location = combined_collection_spec.data["arr"].data_location
        new = combined_collection_spec.replace(arr={"dtype": np.int32})
        assert new.data["arr"].dtype == np.int32
        # Other fields on that DataSetSpec are preserved.
        assert new.data["arr"].data_location == original_location

    def test_replace_data_entry_with_dict_based_on_original_not_copy(
        self, combined_collection_spec
    ):
        """Dict-form data-entry change always applies to the original spec, not the deepcopy."""
        original_dtype = combined_collection_spec.data["arr"].dtype
        new = combined_collection_spec.replace(arr={"dtype": np.int32})
        # The original data entry still reflects the pre-replace dtype.
        assert combined_collection_spec.data["arr"].dtype == original_dtype

    # -----------------------------------------------------------------------
    # Full data-dict replacement (data= keyword)
    # -----------------------------------------------------------------------

    def test_replace_full_data_dict(self, combined_collection_spec, binary_spec):
        """Passing data= as a complete dict replaces the entire data mapping.

        The replacement dict is deep-copied, so the new spec is independent of
        the caller's dict.
        """
        new_data = {"x": binary_spec}
        new = combined_collection_spec.replace(data=new_data)
        assert set(new.data.keys()) == {"x"}
        # Deep-copied: a new DataSetSpec object with the same field values.
        assert new.data["x"] is not binary_spec
        assert new.data["x"].dtype == binary_spec.dtype
        assert new.data["x"].dataset_storage_type == binary_spec.dataset_storage_type
        assert new.data["x"].data_location == binary_spec.data_location

    def test_replace_full_data_dict_combined_with_other_change(
        self, combined_collection_spec, binary_spec
    ):
        """data= replacement can be combined with other field changes."""
        new = combined_collection_spec.replace(
            simulations_combined=False,
            data={"x": binary_spec},
        )
        assert new.simulations_combined is False
        assert set(new.data.keys()) == {"x"}

    # -----------------------------------------------------------------------
    # Multiple simultaneous changes
    # -----------------------------------------------------------------------

    def test_replace_multiple_changes_at_once(
        self, combined_collection_spec, hdf5_attr_spec
    ):
        """Multiple keyword arguments can be applied in a single call."""
        new = combined_collection_spec.replace(
            simulations_combined=False,
            params=hdf5_attr_spec,
            arr={"dtype": np.int32},
        )
        assert new.simulations_combined is False
        assert new.params is hdf5_attr_spec
        assert new.data["arr"].dtype == np.int32

    # -----------------------------------------------------------------------
    # Isolation: original unchanged, result is new frozen instance
    # -----------------------------------------------------------------------

    def test_replace_returns_new_instance(self, combined_collection_spec):
        """replace() always returns a new DataCollectionSpec instance."""
        new = combined_collection_spec.replace(simulations_combined=False)
        assert new is not combined_collection_spec

    def test_replace_does_not_mutate_original_field(self, combined_collection_spec):
        """The original DataCollectionSpec's field is unchanged after replace()."""
        combined_collection_spec.replace(simulations_combined=False)
        assert combined_collection_spec.simulations_combined is True

    def test_replace_does_not_mutate_original_params(self, combined_collection_spec):
        """The original params spec is unchanged after a dict-form params change."""
        original_dtype = combined_collection_spec.params.dtype
        combined_collection_spec.replace(params={"dtype": np.int32})
        assert combined_collection_spec.params.dtype == original_dtype

    def test_replace_result_is_frozen(self, combined_collection_spec):
        """The object returned by replace() is still frozen."""
        new = combined_collection_spec.replace(simulations_combined=False)
        with pytest.raises(FrozenInstanceError):
            new.simulations_combined = True

    # -----------------------------------------------------------------------
    # Unchanged data entries are deep-copied
    # -----------------------------------------------------------------------

    def test_unchanged_data_entries_produce_new_dict(self, combined_collection_spec):
        """When data= is not supplied, the new data dict is a fresh object."""
        new = combined_collection_spec.replace(simulations_combined=False)
        assert new.data is not combined_collection_spec.data

    def test_unchanged_data_entries_are_deep_copied(self, combined_collection_spec):
        """DataSetSpec objects not mentioned in changes are deep-copied."""
        new = combined_collection_spec.replace(simulations_combined=False)
        assert new.data["arr"] is not combined_collection_spec.data["arr"]

    def test_unchanged_data_entry_values_are_preserved(self, combined_collection_spec):
        """Deep-copied unchanged entries keep the same field values."""
        new = combined_collection_spec.replace(simulations_combined=False)
        assert new.data["arr"].dtype == combined_collection_spec.data["arr"].dtype
        assert (
            new.data["arr"].data_location
            == combined_collection_spec.data["arr"].data_location
        )

    # -----------------------------------------------------------------------
    # Error cases
    # -----------------------------------------------------------------------

    def test_unknown_key_raises_key_error(self, combined_collection_spec):
        """A key that is not a field name or a data key raises KeyError."""
        with pytest.raises(KeyError):
            combined_collection_spec.replace(nonexistent_key="value")

    def test_invalid_value_for_field_raises_value_error(
        self, combined_collection_spec
    ):
        """A value that is not a DataSetSpec, bool, or dict raises ValueError."""
        with pytest.raises(ValueError):
            combined_collection_spec.replace(params=42)

    def test_invalid_value_for_data_entry_raises_value_error(
        self, combined_collection_spec
    ):
        """A non-DataSetSpec, non-dict value for a data key raises ValueError."""
        with pytest.raises(ValueError):
            combined_collection_spec.replace(arr="invalid")

    # -----------------------------------------------------------------------
    # Integration: replicates v1.95.0 derivation from v1.99.0
    # -----------------------------------------------------------------------

    def test_replace_on_real_spec_single_field(self):
        """replace() with one dict-form change works on a real dataspec collection."""
        coll = copy.deepcopy(dataspec["v1.99.0"]["microscale_out"])
        new = coll.replace(params={"dtype": np.float64})
        assert new.params.dtype == np.float64
        # Fields not changed are preserved from the original.
        assert new.params.dataset_storage_type == coll.params.dataset_storage_type

    def test_replace_replicates_v1_95_microscale_out_derivation(self):
        """replace() can reproduce the v1.95.0 microscale_out spec from v1.99.0."""
        coll = copy.deepcopy(dataspec["v1.99.0"]["microscale_out"])
        new = coll.replace(
            params={
                "dtype": np.float64,
                "dataset_storage_type": CONST.DATASET_STORAGE_TYPE.FILE_PARSED,
            }
        )
        assert new.params.dtype == np.float64
        assert new.params.dataset_storage_type == CONST.DATASET_STORAGE_TYPE.FILE_PARSED
        # All dataset entries remain intact.
        assert set(new.data.keys()) == set(coll.data.keys())

    def test_replace_result_matches_actual_v1_95_spec(self):
        """Derived microscale_out matches the real v1.95.0 params dtype."""
        v1_95_params = dataspec["v1.95.0"]["microscale_out"].params
        assert v1_95_params.dtype == np.float64
        assert v1_95_params.dataset_storage_type == CONST.DATASET_STORAGE_TYPE.FILE_PARSED
