"""Comprehensive tests for the lysis.data.dataspec module.

Tests cover:
- parse_shape() with integers, variable dims, string references, and errors
- check_dataset_spec() with dtype casting, fixed/variable dims, and params
- DataSetSpec and DataCollectionSpec frozen dataclass behavior
- The dataspec registry dict (versions, tags, collections, dataset names)
- The tags dict (aliases and resolution)
"""

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
