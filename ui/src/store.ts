import { create } from "zustand";

interface PrimerDriverState {
  formKey: number;
  resetForm: () => void;
}

export const usePrimerDriverStore = create<PrimerDriverState>()(set => ({
  formKey: 0,
  resetForm: () => set(state => ({ formKey: state.formKey + 1 })),
}));
