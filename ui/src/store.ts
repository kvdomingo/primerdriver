import { create } from "zustand";

export type Mode = "" | "CHAR" | "DNA" | "PRO";

export type ResultData =
  | Record<string, Record<string, string | number>>
  | { data: string }
  | string;

interface PrimerDriverState {
  loadedResultsFromModule: boolean;
  results: { loaded: boolean; data: ResultData };
  mode: Mode;
  formKey: number;
  setLoadedResultsFromModule: (loaded: boolean) => void;
  setResults: (results: { loaded: boolean; data: ResultData }) => void;
  setMode: (mode: Mode) => void;
  resetForm: () => void;
}

export const usePrimerDriverStore = create<PrimerDriverState>()(set => ({
  loadedResultsFromModule: false,
  results: { loaded: false, data: {} },
  mode: "",
  formKey: 0,
  setLoadedResultsFromModule: loadedResultsFromModule =>
    set({ loadedResultsFromModule }),
  setResults: results => set({ results }),
  setMode: mode => set({ mode }),
  resetForm: () => set(state => ({ formKey: state.formKey + 1 })),
}));
