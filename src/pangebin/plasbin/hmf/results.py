"""HMF result managers."""

from __future__ import annotations

from itertools import product
from typing import TYPE_CHECKING, final

import pangebin.plasbin.bins.input_output as cmn_bins_io
import pangebin.plasbin.bins.items as cmn_bins_items
from pangebin.yaml_interface import YAMLInterface

from . import bins
from . import file_system as fs
from .milp import views as lp_views

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator
    from pathlib import Path


class FeasibleInstance:
    """Feasible bin class instance."""

    @classmethod
    def null(cls) -> FeasibleInstance:
        """Create new feasible instance."""
        return cls(0)

    def __init__(self, number_of_bins: int) -> None:
        self.__number_of_bins = number_of_bins

    def number_of_bins(self) -> int:
        """Get number of bins."""
        return self.__number_of_bins


@final
class BinClass(YAMLInterface):
    """Bin class best feasible instance."""

    KEY_WITH_SEEDS = "with_seeds"
    KEY_FREE_OF_SEEDS = "free_of_seeds"

    @classmethod
    def new(cls) -> BinClass:
        """Create new bin class."""
        return cls(FeasibleInstance.null(), FeasibleInstance.null())

    @classmethod
    def from_dict(cls, obj_dict: dict[str, int]) -> BinClass:
        """Convert dict to object."""
        return cls(
            FeasibleInstance(obj_dict[cls.KEY_WITH_SEEDS]),
            FeasibleInstance(obj_dict[cls.KEY_FREE_OF_SEEDS]),
        )

    def __init__(
        self,
        best_with_seeds: FeasibleInstance,
        best_free_of_seeds: FeasibleInstance,
    ) -> None:
        """Initialize object."""
        self.__best_with_seeds = best_with_seeds
        self.__best_free_of_seeds = best_free_of_seeds

    def best_with_seeds(self) -> FeasibleInstance:
        """Get best feasible instance with seeds."""
        return self.__best_with_seeds

    def best_free_of_seeds(self) -> FeasibleInstance:
        """Get best feasible instance free of seeds."""
        return self.__best_free_of_seeds

    def best_for_seed_constraint(
        self,
        seed_constraint: bins.SeedConstraint,
    ) -> FeasibleInstance:
        """Get best feasible instance from seed constraint."""
        match seed_constraint:
            case bins.SeedConstraint.NOT_REQUIRED:
                return self.best_free_of_seeds()
            case bins.SeedConstraint.REQUIRED:
                return self.best_with_seeds()

    def set_best_with_seeds(self, best_with_seeds: FeasibleInstance) -> None:
        """Set best feasible instance with seeds."""
        self.__best_with_seeds = best_with_seeds

    def set_best_free_of_seeds(self, best_free_of_seeds: FeasibleInstance) -> None:
        """Set best feasible instance free of seeds."""
        self.__best_free_of_seeds = best_free_of_seeds

    def set_best_for_seed_constraint(
        self,
        seed_constraint: bins.SeedConstraint,
        best: FeasibleInstance,
    ) -> None:
        """Set best feasible instance from seed constraint."""
        match seed_constraint:
            case bins.SeedConstraint.NOT_REQUIRED:
                self.set_best_free_of_seeds(best)
            case bins.SeedConstraint.REQUIRED:
                self.set_best_with_seeds(best)

    def to_dict(self) -> dict[str, int]:
        """Convert to dict."""
        return {
            self.KEY_WITH_SEEDS: self.__best_with_seeds.number_of_bins(),
            self.KEY_FREE_OF_SEEDS: self.__best_free_of_seeds.number_of_bins(),
        }


@final
class ConnectedComponent(YAMLInterface):
    """Connected component best feasible instance."""

    KEY_CIRCULAR = "circular"
    KEY_PARTIALLY_CIRCULAR = "partially_circular"

    @classmethod
    def new(cls) -> ConnectedComponent:
        """Create new connected component result."""
        return cls(BinClass.new(), BinClass.new())

    @classmethod
    def from_dict(cls, obj_dict: dict) -> ConnectedComponent:
        """Convert dict to object."""
        return cls(
            BinClass.from_dict(obj_dict[cls.KEY_CIRCULAR]),
            BinClass.from_dict(obj_dict[cls.KEY_PARTIALLY_CIRCULAR]),
        )

    def __init__(
        self,
        best_circular: BinClass,
        best_partially_circular: BinClass,
    ) -> None:
        self.__best_circular = best_circular
        self.__best_partially_circular = best_partially_circular

    def best_circular(self) -> BinClass:
        """Get best feasible instance circular."""
        return self.__best_circular

    def best_partially_circular(self) -> BinClass:
        """Get best feasible instance partially circular."""
        return self.__best_partially_circular

    def best_for_topology(self, topology: bins.Topology) -> BinClass:
        """Get best feasible instance from topology."""
        match topology:
            case bins.Topology.CIRCULAR:
                return self.best_circular()
            case bins.Topology.PARTIALLY_CIRCULAR:
                return self.best_partially_circular()

    def set_best_circular(self, best_circular: BinClass) -> None:
        """Set best feasible instance circular."""
        self.__best_circular = best_circular

    def set_best_partially_circular(self, best_partially_circular: BinClass) -> None:
        """Set best feasible instance partially circular."""
        self.__best_partially_circular = best_partially_circular

    def set_best_for_topology(self, topology: bins.Topology, best: BinClass) -> None:
        """Set best feasible instance for topology."""
        match topology:
            case bins.Topology.CIRCULAR:
                self.set_best_circular(best)
            case bins.Topology.PARTIALLY_CIRCULAR:
                self.set_best_partially_circular(best)

    def to_dict(self) -> dict:
        """Convert to dict."""
        return {
            self.KEY_CIRCULAR: self.__best_circular.to_dict(),
            self.KEY_PARTIALLY_CIRCULAR: self.__best_partially_circular.to_dict(),
        }


@final
class Root(YAMLInterface):
    """Root best results."""

    KEY_CONNECTED_COMPONENTS = "connected_components"

    @classmethod
    def new(cls) -> Root:
        """Create new root result."""
        return cls([])

    @classmethod
    def from_dict(cls, obj_dict: dict[str, list[dict]]) -> Root:
        """Convert dict to object."""
        return cls(
            [
                ConnectedComponent.from_dict(cc)
                for cc in obj_dict[cls.KEY_CONNECTED_COMPONENTS]
            ],
        )

    def __init__(self, connected_components: Iterable[ConnectedComponent]) -> None:
        self.__connected_components = list(connected_components)

    def connected_components(self) -> list[ConnectedComponent]:
        """Get connected components."""
        return self.__connected_components

    def add_connected_component(self, connected_component: ConnectedComponent) -> None:
        """Add connected component."""
        self.__connected_components.append(connected_component)

    def to_dict(self) -> dict[str, list[dict]]:
        """Convert to dict."""
        return {
            self.KEY_CONNECTED_COMPONENTS: [
                cc.to_dict() for cc in self.__connected_components
            ],
        }


# ------------------------------------------------------------------------------------ #
class BinWriter:
    """Result writer for one bin."""

    def __init__(self, bin_fs: fs.Bin) -> None:
        """Initialize object."""
        self.__bin_fs = bin_fs

    def write(
        self,
        bin_stats: bins.Stats,
        milp_stats: lp_views.BinStats,
        seq_normcov: Iterable[cmn_bins_items.SequenceNormCoverage],
    ) -> None:
        """Write bin stats."""
        self.bin_stats(bin_stats)
        self.milp_stats(milp_stats)
        self.seq_normcov(seq_normcov)

    def bin_stats(self, bin_stats: bins.Stats) -> None:
        """Write bin stats."""
        bin_stats.to_yaml(self.__bin_fs.bin_stats_path())

    def milp_stats(self, milp_stats: lp_views.BinStats) -> None:
        """Write MILP stats."""
        milp_stats.to_yaml(self.__bin_fs.milp_stats_path())

    def seq_normcov(
        self,
        seq_normcov: Iterable[cmn_bins_items.SequenceNormCoverage],
    ) -> None:
        """Write bin sequences with their normalized coverages."""
        with cmn_bins_io.Writer.open(self.__bin_fs.bin_seq_normcov_path()) as writer:
            writer.write_bunch_sequences_normcov(seq_normcov)


class BinReader:
    """Result reader for one bin."""

    def __init__(self, bin_fs: fs.Bin) -> None:
        """Initialize object."""
        self.__bin_fs = bin_fs

    def file_system(self) -> fs.Bin:
        """Get bin file system manager."""
        return self.__bin_fs

    def bin_stats(self) -> bins.Stats:
        """Get bin stats."""
        return bins.Stats.from_yaml(self.__bin_fs.bin_stats_path())

    def milp_stats(self) -> lp_views.BinStats:
        """Get MILP stats."""
        return lp_views.BinStats.from_yaml(self.__bin_fs.milp_stats_path())

    def iter_seq_normcov(self) -> Iterator[cmn_bins_items.SequenceNormCoverage]:
        """Get iterator over bin sequences with their normalized coverages."""
        with cmn_bins_io.Reader.open(self.__bin_fs.bin_seq_normcov_path()) as reader:
            yield from reader


class BinClassReader:
    """Bin class result reader."""

    def __init__(
        self,
        feasible_instance: FeasibleInstance,
        bin_class_fs: fs.BinClass,
    ) -> None:
        self.__feasible_instance = feasible_instance
        self.__bin_class_fs = bin_class_fs

    def bin_readers(self) -> Iterable[BinReader]:
        """Get bin readers."""
        for bin_number in range(self.__feasible_instance.number_of_bins()):
            bin_fs = self.__bin_class_fs.bin_file_system(bin_number)
            yield BinReader(bin_fs)

    def file_system(self) -> fs.BinClass:
        """Get bin class file system manager."""
        return self.__bin_class_fs


class ConnectedComponentReader:
    """Connected component result reader."""

    def __init__(
        self,
        best_ccomp_instance: ConnectedComponent,
        ccomp_fs: fs.ConnectedComponent,
    ) -> None:
        self.__best_ccomp_instance = best_ccomp_instance
        self.__ccomp_fs = ccomp_fs

    def circular_with_seeds_reader(self) -> BinClassReader:
        """Get circular bin class reader."""
        return self.reader_from_topology_and_seed_constraint(
            bins.Topology.CIRCULAR,
            bins.SeedConstraint.REQUIRED,
        )

    def circular_free_of_seeds_reader(self) -> BinClassReader:
        """Get circular bin class reader."""
        return self.reader_from_topology_and_seed_constraint(
            bins.Topology.CIRCULAR,
            bins.SeedConstraint.NOT_REQUIRED,
        )

    def partially_circular_with_seeds_reader(self) -> BinClassReader:
        """Get partially circular bin class reader."""
        return self.reader_from_topology_and_seed_constraint(
            bins.Topology.PARTIALLY_CIRCULAR,
            bins.SeedConstraint.REQUIRED,
        )

    def partially_circular_free_of_seeds_reader(self) -> BinClassReader:
        """Get partially circular bin class reader."""
        return self.reader_from_topology_and_seed_constraint(
            bins.Topology.PARTIALLY_CIRCULAR,
            bins.SeedConstraint.NOT_REQUIRED,
        )

    def reader_from_topology_and_seed_constraint(
        self,
        topology: bins.Topology,
        seed_constraint: bins.SeedConstraint,
    ) -> BinClassReader:
        """Get bin class reader from topology and seed constraint."""
        best_feasible_instance = self.__best_ccomp_instance.best_for_topology(
            topology,
        ).best_for_seed_constraint(
            seed_constraint,
        )
        return BinClassReader(
            best_feasible_instance,
            self.__ccomp_fs.bin_class_file_system(
                topology,
                seed_constraint,
                best_feasible_instance.number_of_bins(),
            ),
        )

    def all_bins(self) -> Iterator[BinReaderWithInfo]:
        """Get all bins."""
        return (
            BinReaderWithInfo(
                reader,
                topology,
                seed_constraint,
            )
            for topology, seed_constraint in product(
                (bins.Topology.CIRCULAR, bins.Topology.PARTIALLY_CIRCULAR),
                (bins.SeedConstraint.REQUIRED, bins.SeedConstraint.NOT_REQUIRED),
            )
            for reader in self.reader_from_topology_and_seed_constraint(
                topology,
                seed_constraint,
            ).bin_readers()
        )

    def file_system(self) -> fs.ConnectedComponent:
        """Get connected component file system manager."""
        return self.__ccomp_fs


class RootReader:
    """Root result reader."""

    @classmethod
    def from_output_dir(cls, output_dir: Path) -> RootReader:
        """Get root reader from output directory."""
        root_fs = fs.Root(output_dir)
        best_root_instances = Root.from_yaml(root_fs.best_instances_yaml())
        return cls(best_root_instances, root_fs)

    def __init__(self, best_root_instances: Root, root_fs: fs.Root) -> None:
        """Initialize object."""
        self.__best_root_instances = best_root_instances
        self.__root_fs = root_fs

    def connected_component_readers(self) -> Iterable[ConnectedComponentReader]:
        """Get connected component readers."""
        for ccomp_idx, best_ccomp_instances in enumerate(
            self.__best_root_instances.connected_components(),
        ):
            yield ConnectedComponentReader(
                best_ccomp_instances,
                self.__root_fs.ccomp_file_system(ccomp_idx),
            )

    def all_bins(self) -> Iterator[BinReaderWithInfo]:
        """Get all bins."""
        return (
            bin_rdr_info
            for ccomp_reader in self.connected_component_readers()
            for bin_rdr_info in ccomp_reader.all_bins()
        )


class BinReaderWithInfo:
    """Bin reader with info."""

    def __init__(
        self,
        results_reader: BinReader,
        topology: bins.Topology,
        seed_constraint: bins.SeedConstraint,
    ) -> None:
        """Initialize object."""
        self.__results_reader = results_reader
        self.__topology = topology
        self.__seed_constraint = seed_constraint

    def results_reader(self) -> BinReader:
        """Get bin results reader."""
        return self.__results_reader

    def topology(self) -> bins.Topology:
        """Get topology."""
        return self.__topology

    def seed_constraint(self) -> bins.SeedConstraint:
        """Get seed constraint."""
        return self.__seed_constraint
