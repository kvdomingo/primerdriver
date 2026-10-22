import { useForm, useSelector } from "@tanstack/react-form";
import { useNavigate } from "@tanstack/react-router";
import { MDBCol as Col, MDBRow as Row } from "mdbreact";
import { useEffect } from "react";
import api from "../../api";
import { type ResultData, usePrimerDriverStore } from "../../store";
import { DnaSequenceInput, Form, MutationType, NumberMismatch } from "../form";
import LoadingScreen from "../shared/LoadingScreen";

function Characterize() {
  const navigate = useNavigate();
  const { setMode, setLoadedResultsFromModule, setResults } = usePrimerDriverStore();
  const mode = "CHAR" as const;

  const form = useForm({
    defaultValues: { sequence: "", mismatched_bases: 0, mutation_type: "" },
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
  const isValid = useSelector(
    form.store,
    ({ values }) =>
      values.sequence.length > 0 &&
      values.mismatched_bases > 0 &&
      values.mutation_type !== "",
  );

  useEffect(() => {
    setMode(mode);
  }, [setMode, mode]);

  if (isSubmitting) return <LoadingScreen />;
  return (
    <Form
      title="Primer Characterization"
      onSubmit={form.handleSubmit}
      isValid={isValid}
    >
      <form.Field name="sequence">
        {field => (
          <DnaSequenceInput
            sequence={field.state.value}
            onChange={value =>
              field.handleChange(value.toUpperCase().replace(/[^ACGT]/g, ""))
            }
          />
        )}
      </form.Field>
      <Row className="row-cols-1 row-cols-md-2">
        <Col>
          <form.Field name="mismatched_bases">
            {field => (
              <NumberMismatch
                handleChangeInt={field.handleChange}
                mismatched_bases={field.state.value}
              />
            )}
          </form.Field>
        </Col>
        <Col>
          <form.Field name="mutation_type">
            {field => (
              <MutationType
                handleChange={field.handleChange}
                mutation_type={field.state.value}
              />
            )}
          </form.Field>
        </Col>
      </Row>
    </Form>
  );
}

export default Characterize;
