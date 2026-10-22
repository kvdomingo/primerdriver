import { useForm, useSelector } from "@tanstack/react-form";
import { $api } from "../../api";
import {
  AdvancedSettings,
  DnaSequenceInput,
  Form,
  MutationSelector,
  MutationType,
  defaultSettings,
  designSchema,
} from "../form";
import Result from "./result";

function filterBases(value: string) {
  return value.toUpperCase().replace(/[^ACGT]/g, "");
}

function DnaView() {
  const mode = "DNA" as const;
  const primerDriver = $api.useMutation("post", "/api");

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
    onSubmit: ({ value }) => primerDriver.mutate({ body: { mode, ...value } }),
  });
  const values = useSelector(form.store, state => state.values);
  const isValid = useSelector(
    form.store,
    state => designSchema.safeParse(state.values).success,
  );

  return (
    <>
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
      <Result mode={mode} status={primerDriver.status} data={primerDriver.data} />
    </>
  );
}

export default DnaView;
