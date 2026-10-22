import { useForm, useSelector } from "@tanstack/react-form";
import { useNavigate } from "@tanstack/react-router";
import { useEffect } from "react";
import api from "../../api";
import { type ResultData, usePrimerDriverStore } from "../../store";
import {
  AdvancedSettings,
  DnaSequenceInput,
  Form,
  MutationSelector,
  MutationType,
  defaultSettings,
  designSchema,
} from "../form";
import LoadingScreen from "../shared/LoadingScreen";

function filterBases(value: string) {
  return value.toUpperCase().replace(/[^ACGT]/g, "");
}

function DnaView() {
  const navigate = useNavigate();
  const { setMode, setLoadedResultsFromModule, setResults } = usePrimerDriverStore();
  const mode = "DNA" as const;

  const form = useForm({
    defaultValues: {
      sequence: "",
      mutation_type: "",
      target: "",
      position: 0,
      replacement: "",
      settings: defaultSettings,
    },
    validators: { onSubmit: designSchema },
    onSubmit: async ({ value }) => {
      let data: ResultData;
      try {
        data = (await api.data.primerDriver({ mode, ...value })).data;
      } catch (err) {
        console.error(err);
        data = "Request failed. Please try again later.";
      }
      setLoadedResultsFromModule(true);
      setResults({ data, loaded: true });
      await navigate({ to: "/results" });
    },
  });
  const isSubmitting = useSelector(form.store, state => state.isSubmitting);
  const values = useSelector(form.store, state => state.values);
  const isValid = useSelector(
    form.store,
    state => designSchema.safeParse(state.values).success,
  );

  useEffect(() => {
    setMode(mode);
  }, [setMode, mode]);

  if (isSubmitting) return <LoadingScreen />;
  return (
    <Form
      title="DNA-based Primer Design"
      onSubmit={form.handleSubmit}
      isValid={isValid}
    >
      <form.Field name="sequence">
        {field => (
          <DnaSequenceInput
            sequence={field.state.value}
            onChange={value => field.handleChange(filterBases(value))}
          />
        )}
      </form.Field>
      <form.Field name="mutation_type">
        {field => (
          <MutationType
            handleChange={field.handleChange}
            mutation_type={field.state.value}
          />
        )}
      </form.Field>
      <MutationSelector
        handleChangeTarget={value => form.setFieldValue("target", filterBases(value))}
        handleChangeReplacement={value =>
          form.setFieldValue("replacement", filterBases(value))
        }
        handleChangePosition={value => form.setFieldValue("position", value)}
        mutation_type={values.mutation_type}
        position={values.position}
        replacement={values.replacement}
        target={values.target}
      />
      <form.Field name="settings">
        {field => (
          <AdvancedSettings value={field.state.value} onChange={field.handleChange} />
        )}
      </form.Field>
    </Form>
  );
}

export default DnaView;
